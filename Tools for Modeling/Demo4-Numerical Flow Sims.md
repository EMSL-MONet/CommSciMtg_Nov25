# Tools for Modeling
## Demo 4: Performing Numerical Simulations on Example Soil Micromodel
## Workshop Session
What: Session 4: In-Person Training   <br>
When: Wednesday, February 5th: 1:00 pm – 2:30 am, or  2:40 pm - 4:10 pm <br>
Where:  EMSL 1077, Breakout Room 2

**Lal Mamud** | Postdoctoral Researcher   <br>

## Know Before You Go
This in-person training session is intended to provide an over tools. The material for each of these demonstrations is provided within folders and files of this repo. Participants should make every effort to utilize the scripts provided in order to get the most out of the session, even if they are more familiar with another scripting languages.
### Pre-session Checklist:
- Test Your Environment: Test a few code blocks in the Google Colab link provided below to ensure everything is running smoothly.
- Repository Overview: Familiarize yourself with the materials in the repository to understand the workflow and tools being used.

## Data Access
All the materials needed for the hands-on session are located in this repository: https://github.com/EMSL-MONet/CommSciMtg_Nov25/tree/main/Tools%20for%20Modeling

### Learning Objective
After completing this in-person training, you will be able to:
- Execute and interpret workflows using Jupyter Notebooks on Google Colab.
- Mapping permeability field for the micromodels of interst.
- Understand the fundamentals of numerical simulations applied to soil micromodels.
- Simulate 2D steady-state fulid flow in soil micromodels.

### Before The Day
- Check the readmes beforehand!
- Open the Notebook on Google Colab 🌐 – Log in with your Google account and set up early to avoid delays.
- Test It 🛠️ – Run the first few code blocks to ensure everything’s working smoothly. If it doesn't run let one of us know. 

### On The Day
- Bring Your Laptop 💻 – Ensure it is connected to the internet and ready for use.
- Arrive a few minutes early to get settled and ensure your environment is functioning as expected.

## Agenda on contents of the workshop
Use Jupyter Notebook deployed on Google Colab to execute the workflow.  
  - **Link to collab JUPYTER NOTEBOOK DEMO 4 https://colab.research.google.com/drive/1RJMzxUXDq_602R1egBkcGNTSfXNHdPIN?usp=sharing**
  
  - **Link to toy data for notebooks https://github.com/aramyxt/MONet_Pore2Chip_data**

This Jupyter Notebook may require a high amount of RAM to run depending on the mesh size. The default (free) plan of Google Colab may not provide sufficient memory. 
If you encounter memory issues, please follow these steps to run the notebook locally:
- Set Up a Local Environment: Create a Conda or Anaconda environment to ensure compatibility and easy package management.
- Install Required Packages:
  - Install the pore2chip using
      ```
      pip install pore2chip
      ```
  - Install Jupyter Notebook to run the script locally using
      ```
      pip install notebook
      ```
- Download the Notebook and Data: Download the notebook (p2c_example_2_flow_2d_numerical_on_chip.ipynb) and associated data using the provided links.
- Run Locally: Open the downloaded notebook in Jupyter Notebook from your local directory and execute the cells.
