# Instructions for using WiMBCode
1. Run startup.m.
2. Create a folder containing model equations and parameters to use with a dde solver. The spreadsheet containing parameter values should list parameter names in the first column and values in the second column. See wright20/wright_pars.xlsx for an example.
3. Add datasets to the data folder following the same format as in data/WeltData.xlsx. The program currently imports sheets  labeled 'LH', 'FSH', 'E2', 'P4'.
4. Add new folder(s) to startup.m.
5. Run simulations using the readrun function with inputs specific to the model. Alternatively, use main.m to plot solutions.
