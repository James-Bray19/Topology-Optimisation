# Topology Optimisation with FEA

## Introduction

This application is designed to help you perform topology optimization on 2D models. In this readme, you'll find essential information on how to use the app effectively.

![270122948-95a0b44f-8cf2-440d-88c1-7b03331ccf8e](https://github.com/James-Bray19/Topology-Optimisation/assets/47334864/187feda0-55b9-411f-8ac5-75af1fcfdc84)

## Features

### Topology Optimization

The primary purpose of this application is to perform topology optimization. Topology optimization is a technique that determines the optimal distribution of material within a given design space, subject to certain constraints. In this app:

- **Constraint Editor**: You can set up constraints by moving and editing them using the options panel on the right side of the screen.

- **Volume Fraction**: Define the desired volume fraction, which represents the amount of material you want to allocate within the model.

- **Optimization**: Once constraints and volume fraction are set, simply hit the "Optimise Topology" button to initiate the optimization process.

- **Optimal Topology**: The resulting optimal topology will be displayed in the bottom right corner of the screen. You can choose your preferred view mode (mesh or pixel) and coloring (color or grayscale) from the options on the right.

### Constraint Presets

- **Save and Load Presets**: You can save your constraint configurations and load them later using the buttons provided in the toolbar. This feature allows you to reuse your setups efficiently.

### Displacement Visualizer

- **Displacement Visualization**: The app also provides a displacement visualizer. It renders the model as a mesh and shows how the model moves under specific conditions. You can control the amount of displacement through the options panel on the right.

## Getting Started

To get started with the Unity 2D Topology Optimization App, follow these steps:

1. Download the files and launch this file: `Builds\Windows\Version 1\Topology Optimisation with FEA.exe`

2. Set up your constraints using the Constraint Editor on the top left of the screen.

3. Define the desired volume fraction.

4. Click the "Optimise Topology" button to initiate the optimization process.

5. Observe the optimal topology in the bottom right corner of the screen.

6. Customize the view mode and color as needed.

7. Optionally, save your constraint configurations for future use.

8. Explore the displacement visualizer to understand how the model behaves under different conditions.

## Important Notes

- **Single Iteration Optimization**: Please note that this software performs only one iteration of the optimization process due to time constraints, as it was developed as part of an A-level computer science project. If you require more extensive optimization, you may need to explore other tools or software, or improve mine and submit a pull request.

- **Help Button**: If you're unsure how to use the application, there is a "Help" button located in the top-right corner. This feature provides guidance on using most of the software's functionalities.

- **Code Examples**: These code samples are scripts designed to work inside a unity environment and will not run on their own.

## Feedback and Support

If you encounter any issues, have suggestions for improvements, or need assistance, please feel free to reach out to me.

To edit the Unity project, please contact me for the files, since they were too large for the repo.
