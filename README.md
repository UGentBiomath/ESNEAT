
# Enhancing Accuracy and Efficiency in Calibration of Drinking Water Distribution Networks Through Evolutionary Artificial Neural Networks and Expert Systems

This repository, **ESNEAT**, contains the implementation of the automatic (re)calibration method **ES-NEAT** for the Drinking water distribution network of Fossolo. The approach is described in the scientific paper *"Enhancing Accuracy and Efficiency in Calibration of Drinking Water Distribution Networks Through Evolutionary Artificial Neural Networks and Expert Systems"*. This project is targeted at researchers interested in model calibration, NEAT algorithms, and water resource management.

---

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Input Data](#input-data)
- [Results](#results)
- [Citing This Work](#citing-this-work)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [Contact Information](#contact-information)

---

## Overview
This project aims to enhance reproducibility and transparency in dwdn modeling by implementing **ES-NEAT** to automatically (re)calibrate model parameters, ensuring continuous model performance. The DWDN model and data used in this project are based on the Fossolo network.

### Key Concepts
- **NEAT (NeuroEvolution of Augmenting Topologies)**: A genetic algorithm for evolving neural networks.
- **Expert Systems**: Used to guide the recalibration process.
- **Calibration Methods**: Ensuring model outputs align with real-world data.

---

## Features
- Automatic adjustment of DWDN parameters.
- Visualization of results comparing simulation outputs to real-world data.
- Model reproducibility with included NEAT seed files.
- Python implementation with minimal dependencies.

---

## Installation
### Dependencies
To run this project, you will need the following:
- Python 3.12
- Numpy 1.26.4
- Pandas 2.2.2
- NEAT 0.92
- WNTR 1.4.0
- EPANET 2.2 (required for additional DWDN modeling functionality)

### Setup
1. Clone this repository:
   ```bash
   git clone https://github.com/your-repo-url/ESNEAT.git
   cd ESNEAT
   ```
---

## Usage
### Tools
You can run the code using:
- **Spyder**
- **Jupyter Notebook**

### Steps to Run
1. Load the input data from the `input_Data` folder.
2. Execute the calibration script or notebook.
3. View the generated results including RMSE, parameter adjustments, and visual plots.

No special environment configurations are required.

---

## Input Data
The repository includes the necessary input data located in the `input_Data` folder. The data is provided in:
- **CSV files** for tabular inputs.
- **TXT files** for additional configuration details.

---

## Results
The code generates the following outputs:
- Root Mean Square Error (RMSE) of effluent data compared to real-world measurements.
- Visual plots of flow and pressure data.
- Final adjusted model parameters.
- Time taken for the simulation.
- Updated NEAT seed file for future recalibrations.

Example results are detailed in the accompanying paper.

---

## Citing This Work
*Citation information will be added here when available.*

---

## Contributing
Contributions are welcome! If you'd like to contribute, please send your suggestions or improvements via email.

---

## License
This project is licensed under the **MIT License**. Refer to the `LICENSE` file for more details on usage restrictions and permissions.

---

## Acknowledgments
*Acknowledgment details will be added in the future.*

---

## Contact Information
For questions or support, please refer to the contact information provided in the paper.
