using System.Collections;
using System.Collections.Generic;
using static System.Math;
using UnityEngine;
using UnityEngine.UI;
using MathNet.Numerics;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Factorization;

public class MarchingSquares : MonoBehaviour
{
    #region Variables
    public int xSize;
    public int ySize;
    public Vector<double> U;
    public Gradient gradient;
    public Slider dispExg;

    [HideInInspector] // holds the indices of vertices for each triangle
    private static List<int[]> triTable = new List<int[]>() { 
        new int[] { },                           //0000
        new int[] { 7, 6, 3 },                   //0001
        new int[] { 6, 5, 2 },                   //0010
        new int[] { 3, 7, 5, 3, 5, 2 },          //0011
        new int[] { 7, 0, 4 },                   //0100
        new int[] { 3, 0, 4, 3, 4, 6 },          //0101
        new int[] { 7, 0, 4, 2, 6, 5 },          //0110
        new int[] { 3, 0, 4, 3, 4, 5, 3, 5, 2 }, //0111
        new int[] { 5, 4, 1 },                   //1000
        new int[] { 7, 6, 3, 5, 4, 1 },          //1001
        new int[] { 6, 4, 2, 2, 4, 1 },          //1010
        new int[] { 3, 7, 2, 2, 7, 4, 2, 4, 1 }, //1011
        new int[] { 7, 0, 1, 7, 1, 5 },          //1100
        new int[] { 0, 1, 5, 0, 5, 6, 0, 6, 3 }, //1101
        new int[] { 1, 2, 6, 1, 6, 7, 1, 7, 0 }, //1110
        new int[] { 0, 1, 3, 1, 2, 3 } };        //1111
    #endregion

    public void Render(Matrix<double> p, bool isColored)
    {
        // if theres a displacement vector passed in, interpolate it between 0 and 1
        if (U != null) 
        {
            double max = 0;
            for (int i = 0; i < U.Count; i++) { if (Abs(U[i]) > max) { max = Abs(U[i]); } }
            for (int i = 0; i < U.Count; i++) { U[i] /= max; }
            U *= dispExg.value;
        }

        // combines every element mesh
        CombineInstance[] combine = new CombineInstance[(xSize) * (ySize)];
        int meshIndex = 0;
        for (int i = 0; i < xSize; i++)
        {
            for (int j = 0; j < ySize; j++)
            {
                combine[meshIndex].mesh = GetElementMesh(i, j, p);
                combine[meshIndex].transform = Matrix4x4.Translate(new Vector3(i + 0.5f, j + 0.5f, 0f));
                meshIndex++;
            }
        }

        // renders mesh
        GetComponent<MeshFilter>().mesh = new Mesh();
        GetComponent<MeshFilter>().mesh.CombineMeshes(combine);

        // sets the corner of each mesh to a color based on stress
        // a vertex shader handles the color gradient
        if (U == null && isColored)
        {
            Mesh chunkMesh = GetComponent<MeshFilter>().mesh;
            Color[] colors;
            colors = new Color[chunkMesh.vertices.Length];
            for (int i = 0; i < colors.Length; i++)
            {
                colors[i] = gradient.Evaluate((float)p[(int)chunkMesh.vertices[i].x, (int)chunkMesh.vertices[i].y]);
            }
            chunkMesh.colors = colors;
        }
    }

    public Mesh GetElementMesh(int x, int y, Matrix<double> p)
    {
        bool[] arr = new bool[4];

        //for each vertex
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                // generate four bits, one for each vertex
                if (p[x+i, y+j] >= 0.5f) { arr[i+2*j] = true; } // true = inside
                else { arr[i+2*j] = false; }
            }
        }

        Vector3[] relativePoints = new Vector3[8] 
        { new Vector3 (-0.5f, 0.5f, 0f), new Vector3 ( 0.5f, 0.5f, 0f), // vertices
          new Vector3 ( 0.5f,-0.5f, 0f), new Vector3 (-0.5f,-0.5f, 0f), // 
          new Vector3 (0f, 0.5f, 0f),    new Vector3 ( 0.5f, 0f, 0f),   // midpoints
          new Vector3 (0f,-0.5f, 0f),    new Vector3 (-0.5f, 0f, 0f)};  //

        //if the element is on the edge of the mesh, interpolate midpoints
        if ((arr[0] || arr[1] || arr[2] || arr[3]) && !(arr[0] && arr[1] && arr[2] && arr[3]))
        {
            relativePoints[4] = new Vector3(Mathf.InverseLerp((float)p[x, y + 1], (float)p[x + 1, y + 1], 0.5f) - 0.5f, 0.5f, 0f);
            relativePoints[5] = new Vector3(0.5f, Mathf.InverseLerp((float)p[x + 1, y], (float)p[x + 1, y + 1], 0.5f) - 0.5f, 0f);
            relativePoints[6] = new Vector3(Mathf.InverseLerp((float)p[x, y], (float)p[x + 1, y], 0.5f) - 0.5f, -0.5f, 0f);
            relativePoints[7] = new Vector3(-0.5f, Mathf.InverseLerp((float)p[x, y], (float)p[x, y + 1], 0.5f) - 0.5f, 0f);
        }

        if (U != null)
        {
            List<float> u = eDisps(x, y);
            relativePoints[0] += new Vector3(u[0], u[1], 0f);
            relativePoints[1] += new Vector3(u[2], u[3], 0f);
            relativePoints[2] += new Vector3(u[4], u[5], 0f);
            relativePoints[3] += new Vector3(u[6], u[7], 0f);

            //interpolate
            relativePoints[4] += new Vector3(((u[0] + u[2]) / 2), ((u[1] + u[3]) / 2), 0f);
            relativePoints[5] += new Vector3(((u[2] + u[4]) / 2), ((u[3] + u[5]) / 2), 0f);
            relativePoints[6] += new Vector3(((u[4] + u[6]) / 2), ((u[5] + u[7]) / 2), 0f);
            relativePoints[7] += new Vector3(((u[6] + u[0]) / 2), ((u[7] + u[1]) / 2), 0f);
        }

        // convert to integer
        BitArray bitArray = new BitArray(arr);
        int[] array = new int[1];
        bitArray.CopyTo(array, 0);

        // output the pixel mesh
        Mesh voxelMesh = new Mesh();
        voxelMesh.vertices = relativePoints;
        voxelMesh.triangles = triTable[array[0]];

        if (U != null)
        {
            Color[] colors;
            colors = new Color[voxelMesh.vertices.Length];
            for (int i = 0; i < colors.Length; i++)
            {
                if (((x + y) % 2) == 0) { colors[i] = new Color(0.6509434f, 0.4485452f, 0f, 1f); }
                else { colors[i] = new Color(0.4528302f, 0.3109798f, 0f, 1f); }
            }
            voxelMesh.colors = colors;
        }

        return voxelMesh;
    }

    List<float> eDisps(int x, int y)
    {
        int n1 = 2 * ((ySize + 2) * x + 1 + ySize - y);
        int n2 = 2 * ((ySize + 2) * (x + 1) + 1 + ySize - y);
        int[] e = new int[] { n1 - 2, n1 - 1, n2 - 2, n2 - 1, n2, n2 + 1, n1, n1 + 1 };

        List<float> l = new List<float>();
        foreach (int i in e) { l.Add((float)U[i]); }
        return l;
    }
}