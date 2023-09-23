using System.Collections;
using System.Collections.Generic;
using System.Linq;
using static System.Math;
using UnityEngine;
using Unity.Collections;
using UnityEngine.UI;
using MathNet.Numerics;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Factorization;

public class Optimiser : MonoBehaviour
{
    #region Attributes
    public int xSize, ySize; // domain size
    public static string optimiseStatus; // optimisation stage
    public int iterations; // temp
    private int iteration; // temp

    public double penal; // penalisation factor
    public double rMin; // density filter radius
    public double pMin; // smallest density value
    public double tolerance; // smallest density value
    public double move; // move limit

    private Matrix<double> LSM; // relates element force to displacement
    private Matrix<double> GSM; // relates force to displacement
    private Matrix<double> p; // densities calculated from potential energy

    private Vector<double> U; // displacement vector
    private Vector<double> F; // force vector


    public RawImage topoTexture; // output image for texture
    public Toggle topoMesh, topoColor; // post processing controls
    public MarchingSquares topoMarcher; // marching cubes script
    public MarchingSquares dispMarcher; // marching cubes script

    public Slider volumeSlider; // user input slider for the volume fraction
    public Gradient gradient; // colour gradient to display with
    public GraphGenerator graph; // graph script to send values to

    public GameObject backToEdit;
    #endregion


    // encapsulation so UI components can change the status
    public void SetOptimiseStatus(string status)
    {
        optimiseStatus = status;
    }


    // called at the start of runtime
    void Start()
    {
        optimiseStatus = "idle";
        LSM = DenseMatrix.OfArray(new double[,] {
        { 0.4945054945054945,0.1785714285714285,-0.3021978021978022,-0.01373626373626375,-0.2472527472527473,-0.1785714285714285,0.05494505494505494,0.01373626373626375 },
        { 0.1785714285714285,0.4945054945054945,0.01373626373626375,0.05494505494505494,-0.1785714285714285,-0.2472527472527473,-0.01373626373626375,-0.3021978021978022 },
        { -0.3021978021978022, 0.01373626373626375, 0.4945054945054945, -0.1785714285714285, 0.05494505494505494, -0.01373626373626375, -0.2472527472527473, 0.1785714285714285},
        { -0.01373626373626375,0.05494505494505494,-0.1785714285714285,0.4945054945054945,0.01373626373626375,-0.3021978021978022,0.1785714285714285,-0.2472527472527473},
        { -0.2472527472527473,-0.1785714285714285,0.05494505494505494,0.01373626373626375,0.4945054945054945,0.1785714285714285,-0.3021978021978022,-0.01373626373626375},
        { -0.1785714285714285,-0.2472527472527473,-0.01373626373626375,-0.3021978021978022,0.1785714285714285,0.4945054945054945,0.01373626373626375,0.05494505494505494},
        { 0.05494505494505494,-0.01373626373626375,-0.2472527472527473,0.1785714285714285,-0.3021978021978022,0.01373626373626375,0.4945054945054945,-0.1785714285714285},
        { 0.01373626373626375,-0.3021978021978022,0.1785714285714285,-0.2472527472527473,-0.01373626373626375,0.05494505494505494,-0.1785714285714285,0.4945054945054945}});
    }


    // called every frame
    void Update()
    {

        if (optimiseStatus == "iterating")
        {
            if (iteration > 0)
            {
                graph.Clear(); // temp
                F = GetForceVector(); // assemble force vector
                GSM = GetGlobalMatrix(); // assemble GSM
                ReduceEquation(); // remove fixpoints
                U = GSM.Solve(F); // solve linear system
                ExpandResult(); // add fixpoints

                p = GetPotentials(); // set densities
                CriteriaUpdate(); // sensitivity analysis, filter, then match criteria

                Display();
                iteration--;
            }
            else { optimiseStatus = "optimised"; }
        }
        if (optimiseStatus == "optimised") { backToEdit.SetActive(true); }
    }


    #region FEA
    private Vector<double> GetForceVector()
    {
        Vector<double> f = new DenseVector((xSize + 1) * (ySize + 1) * 2);
        foreach (GameObject force in InputManager.forces)
        {
            Vector2 pos = force.GetComponent<RectTransform>().anchoredPosition;
            int n = ndof((int)pos.x, (int)pos.y);
            f[n] = force.GetComponentInChildren<ForceManager>().force.x;
            f[n + 1] = force.GetComponentInChildren<ForceManager>().force.y;
        }
        return f;
    }


    // assembles the global stiffness matrix
    private Matrix<double> GetGlobalMatrix()
    {
        int matrixSize = (xSize + 1) * (ySize + 1) * 2;
        Matrix<double> K = new DenseMatrix(matrixSize, matrixSize);

        // for each element
        for (int j = 0; j < ySize; j++)
        {
            for (int i = 0; i < xSize; i++)
            {
                int[] edof = edofs(i, j);
                int xInd = 0, yInd = 0;
                foreach (int y in edof)
                {
                    foreach (int x in edof)
                    {
                        // add the DOF into the GSM weighted by density
                        K[x, y] += p[i, j] * LSM[xInd, yInd];
                        xInd++;
                    }
                    yInd++;
                    xInd = 0;
                }
            }
        }
        return K;
    }


    // remove fixed DOFs from force vector and GSM
    private void ReduceEquation()
    {
        // convert to matrix as vector doesn't support row/column removal
        Matrix<double> F_asColMat = F.ToColumnMatrix();
        List<int> toRemove = GetFixpointIndices();
        for (int i = toRemove.Count - 1; i >= 0; i--)
        {
            F_asColMat = F_asColMat.RemoveRow(toRemove[i]);
            GSM = GSM.RemoveRow(toRemove[i]).RemoveColumn(toRemove[i]);
        }
        F = new DenseVector(F_asColMat.ToColumnArrays()[0]);
    }


    // add zeros in place of the fixpoints
    private void ExpandResult()
    {
        Vector<double> zero = new DenseVector(1);
        List<double> U_asList = U.ToArray().ToList();
        List<int> toAdd = GetFixpointIndices();
        for (int i = 0; i < toAdd.Count; i++)
        {
            if (toAdd[i] >= U_asList.Count) {  U_asList.Add(0); }
            else { U_asList.Insert(toAdd[i], 0); }
        }
        U = new DenseVector(U_asList.ToArray());
    }


    // using the displacement vector, compute elemental energy, then output density
    Matrix<double> GetPotentials()
    {
        Matrix<double> temp = new DenseMatrix(xSize, ySize);
        Vector<double> u = new DenseVector(8); // holds the displacements of the elements

        for (int j = 0; j < ySize; j++)
        {
            for (int i = 0; i < xSize; i++)
            {
                // adds the displacements to the element displacement vector
                int[] edof = edofs(i, j);
                for (int x = 0; x < 8; x++) { u[x] = U[edof[x]]; }

                // calculate stress
                temp[i, j] = (u.ToRowMatrix() * LSM * u)[0];
            }
        }
        return temp;
    }
    #endregion


    #region Optimisation
    // starts process
    public void Optimise() 
    {
        if (ConditionsAreValid())
        {
            foreach (GameObject f in InputManager.forces) { f.GetComponentInChildren<Draggable>().enabled = false; }
            foreach (GameObject f in InputManager.fixpoints) { f.GetComponentInChildren<Draggable>().enabled = false; }

            iteration = iterations;
            graph.Clear();
            InputManager.activeObject = null;

            p = Matrix<double>.Build.Dense(xSize, ySize, volumeSlider.value);
            optimiseStatus = "iterating";
        }
        else 
        {
            optimiseStatus = "optimised";
        }
    }


    // apply a 'blur' type effect to eliminate checkerboarding and anomalies
    Matrix<double> Filter(Matrix<double> matrix)
    {
        Matrix<double> temp = new DenseMatrix(xSize, ySize);
        for (int j = 0; j < ySize; j++)
        {
            for (int i = 0; i < xSize; i++)
            {
                List<double> weights = new List<double>();
                for (int v = -(int)rMin; v < (int)rMin; v++)
                {
                    if (j + v > 0 && j + v < ySize)
                    {
                        for (int u = -(int)rMin; u < (int)rMin; u++)
                        {
                            if (i + u > 0 && i + u < xSize)
                            {
                                double dist = Sqrt(v * v + u * u);
                                if (dist < rMin)
                                {
                                    //weights.Add((1 - (dist / rMin)) * matrix[i + u, j + v]);
                                    weights.Add((1 - (dist / rMin)) * matrix[i + u, j + v]);
                                }
                            }
                        }
                    }
                }
                temp[i, j] = weights.Average();
            }
        }
        return temp;
    }

    //Matrix<double> Filter(Matrix<double> matrix)
    //{
    //    Matrix<double> temp = new DenseMatrix(xSize, ySize);
    //    for (int i = 0; i < xSize; i++)
    //    {
    //        for (int j = 0; j < ySize; j++)
    //        {
    //            double sum = 0;
    //            for (int u = -(int)rMin; u < (int)rMin; u++)
    //            {
    //                if (i + u >= 0 && i + u < xSize)
    //                {
    //                    for (int v = -(int)rMin; v < (int)rMin; v++)
    //                    {
    //                        if (j + v >= 0 && j + v < ySize)
    //                        {
    //                            double fac = Max(0, rMin - Sqrt(v * v + u * u));
    //                            sum = sum + fac;
    //                            temp[i, j] = temp[i, j] + fac * p[u, v] * matrix[u, v];
    //                        }
    //                    }
    //                }
    //            }
    //            temp[i, j] = temp[i, j] / (p[i, j] * sum);
    //        }
    //    }
    //    return temp;
    //}


    // bisection algorithm to meet the users volume criteria
    void CriteriaUpdate()
    {
        // get dc matrix
        Matrix<double> dc = new DenseMatrix(xSize, ySize);
        for (int j = 0; j < ySize; j++) 
        {
            for (int i = 0; i < xSize; i++) 
            {
                dc[i, j] = -penal * Pow(p[i, j], penal - 1) * p[i, j];
            }
        }

        // filter dc matrix
        dc = Filter(dc);

        // update to match criteria
        Matrix<double> pTrial = p;
        double lambda, lower = 0, upper = 2.3424;
        for (int attempts = 0; attempts < 100; attempts++)
        {
            if (Abs(pTrial.RowSums().Sum() / (xSize * ySize) - volumeSlider.value) < tolerance) { break; }

            lambda = (upper + lower) / 2;
            for (int j = 0; j < ySize; j++)
            {
                for (int i = 0; i < xSize; i++)
                {
                    pTrial[i, j] = Max(pMin, Max(p[i, j] - move, Min(1f, Min(p[i, j] + move, p[i, j] * Sqrt(-dc[i, j] / lambda)))));
                }
            }
            //graph.AddPoint(pTrial.RowSums().Sum() / (xSize * ySize)); // temp
            if (pTrial.RowSums().Sum() / (xSize * ySize) > volumeSlider.value) { lower = lambda; } else { upper = lambda; }
        }
        p = pTrial;
    }


    // return a matrix where every value is mapped between 0 and 1
    //Matrix<double> Interpolate(Matrix<double> matrix)
    //{
    //    double highestValue = 0;
    //    for (int j = 0; j < matrix.ColumnCount; j++)
    //    {
    //        for (int i = 0; i < matrix.RowCount; i++)
    //        {
    //            if (matrix[i, j] > highestValue) { highestValue = matrix[i, j]; }
    //        }
    //    }
    //    for (int j = 0; j < matrix.ColumnCount; j++)
    //    {
    //        for (int i = 0; i < matrix.RowCount; i++)
    //        {
    //            double val = matrix[i, j];
    //            val = ((matrix[i, j] * (1f - pMin)) / highestValue) + pMin;
    //            matrix[i, j] = val;
    //        }
    //    }
    //    return matrix;
    //}
    #endregion


    #region Index Functions
    // returns the DOF indices for an element
    int[] edofs(int x, int y)
    {
        int n1 = ndof(x, y);
        int n2 = ndof(x + 1, y);
        //string output = "";
        //foreach (int k in (new int[] { n1 - 2, n1 - 1, n2 - 2, n2 - 1, n2, n2 + 1, n1, n1 + 1 })) { output += k.ToString() + ", "; }
        //Debug.Log("edof for element x: " + x.ToString() + " y: " + y.ToString() + " is " + output);
        return new int[] { n1 - 2, n1 - 1, n2 - 2, n2 - 1, n2, n2 + 1, n1, n1 + 1 };
    }


    // returns the DOF index of a node
    int ndof(int x, int y)
    {
        return 2 * ((ySize + 1) * x + ySize - y);
    }


    // returns a sorted list of DOFs for all fixpoints
    List<int> GetFixpointIndices()
    {
        List<int> indices = new List<int>();
        foreach (GameObject f in InputManager.fixpoints)
        {
            Vector2 pos = f.GetComponent<RectTransform>().anchoredPosition;
            if (f.GetComponentInChildren<FixpointManager>().x) { indices.Add(ndof((int)pos.x, (int)pos.y)); }
            if (f.GetComponentInChildren<FixpointManager>().y) { indices.Add(ndof((int)pos.x, (int)pos.y) + 1); }
        }
        indices.Sort();
        //string output = "";
        //foreach(int i in indices) { output += " " + i; }
        //Debug.Log(output);
        return indices;
    }
    #endregion


    // checks there are no errors with the inputs
    bool ConditionsAreValid()
    {
        foreach (GameObject f in InputManager.forces)
        {
            f.GetComponentInChildren<Draggable>().enabled = false;
            if (!f.GetComponentInChildren<ForceManager>().isValid) 
            {
                ErrorText.ShowError("Invalid force detected");
                return false; 
            }
        }
        foreach (GameObject f in InputManager.fixpoints)
        {
            f.GetComponentInChildren<Draggable>().enabled = false;
            if (!f.GetComponentInChildren<FixpointManager>().isValid) 
            {
                ErrorText.ShowError("Invalid fixpoint detected");
                return false; 
            }
        }

        if (InputManager.forces.Count == 0) 
        {
            ErrorText.ShowError("Force required");
            return false; 
        }
        else if (InputManager.fixpoints.Count < 2) 
        {
            ErrorText.ShowError("2 or more fixpoints required");
            return false; 
        }
        return true;
    }


    #region Outputs
    // output model defined by post processing toggles
    public void Display()
    {
        dispMarcher.GetComponent<MeshRenderer>().enabled = true;
        dispMarcher.U = U;
        dispMarcher.Render(p, false);

        if (topoMesh.isOn)
        {
            topoTexture.enabled = false;
            topoMarcher.GetComponent<MeshRenderer>().enabled = true;
            topoMarcher.U = null;
            topoMarcher.Render(p, topoColor.isOn);
        }
        else
        {
            topoTexture.enabled = true;
            topoMarcher.GetComponent<MeshRenderer>().enabled = false;
            Texture2D texture = new Texture2D(xSize, ySize);
            topoTexture.texture = texture;
            for (int j = 0; j < ySize; j++)
            {
                for (int i = 0; i < xSize; i++)
                {
                    if (topoColor.isOn) { texture.SetPixel(i, j, gradient.Evaluate((float)p[i, j])); }
                    else { texture.SetPixel(i, j, new Color((float)p[i, j], (float)p[i, j], (float)p[i, j], 1f)); }
                }
            }
            texture.filterMode = FilterMode.Point;
            texture.Apply();
        }
    }

    public void UpdateDisplay()
    {
        if (optimiseStatus != "idle") { Display(); }
    }
    #endregion
}