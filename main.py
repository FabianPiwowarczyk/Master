import netCDF4 as nc
import IASI


file_path = "data/IASIA_MUSICA_030300_L2_AllTargetProducts_20200101000557_68496.nc"

# Open the dataset
dataset = nc.Dataset(file_path, mode='r')

# Print a list of all variable names
print("Variable names:")
for var_name in dataset.variables:
    print(var_name)


def printvar(var_name):
    # Access the variable
    variable = dataset.variables[var_name]

    # Print the variable's name
    print("Variable:", var_name)

    # Print the variable's description (if available)
    if hasattr(variable, 'description'):
        print("Description:", variable.description)

    # Print the variable's data
    print("Data:")
    print(variable[:])  # Print all the values of the variable

    # If the variable has dimensions, you might want to print them as well
    if hasattr(variable, 'dimensions'):
        print("Dimensions:", variable.dimensions)


printvar('musica_wv')

#print(IASI.read_all_iasi('data'))
