
# About the Project

This repository features MATLAB code designed for constructing lineage trees from time-lapse imaging of bacterial cells. The trees represent cell division history and annotate cells with various characteristics, such as cell type based on fluorescence intensity. This is achieved by analyzing data typically obtained from cell tracking and segmentation software like CellProfiler, with each cell's metrics stored in separate Excel files for different microfluidic device chambers.

This project significantly relies on the foundational work provided by [Jean-Yves Tinevez](http://tinevez.github.io/matlab-tree/) for tree construction algorithms and data structures. Our code extends this work to apply specifically to bacterial cell lineage tracking. Proper acknowledgment and appreciation are extended to this original work, which serves as a critical component of our project's functionality.

To successfully run this code:

1. You must download and integrate Tinevez's MATLAB tree construction toolkit available at [this repository](http://tinevez.github.io/matlab-tree/). Ensure it is downloaded and placed within the same working directory as this project, or appropriately referenced within your MATLAB path settings.

2. Prepare your segmented and tracked cell data in Excel format, which should be an output from CellProfiler or similar segmentation software.

3. The main script, `Main.m`, is the entry point for the workflow. It requires customization to align with your dataset, particularly linking to your Excel files and ensuring data headers match your file's format.

4. Key modifications and data path settings should be adjusted within `Main.m` and the `BuildTree.m` class to accommodate your specific dataset and tracking information.


# Contributions

I welcome feedback, suggestions, and contributions to enhance this project. Should you encounter issues or wish to propose improvements, feel free to open an issue or submit a pull request.

### Contact

For further information, questions, or collaboration proposals, please don't hesitate to get in touch via email at sshaghay@uwaterloo.ca .
