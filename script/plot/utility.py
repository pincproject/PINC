import numpy as np

# Takes a hf-data and transforms it into 1D - array
# retains x axis.
def transformData(dim,dataset, timestep):
    grid = dataset['/n=%.1f'%timestep]
    grid = np.squeeze(grid)
    if len(grid.shape)==4:
        grid = grid[:,:,:,0]
    grid = np.transpose(grid) #Accounting for reverse zyx order
    nDims = len(grid.shape)
    for d in range(nDims):
        if d!=dim:
            if d<dim:
                grid = np.average(grid, axis = 0)
            else:
                grid = np.average(grid, axis = 1)


    return grid

#   Plots a 1D subgrid
#   input:  string      name
#           1D array    grid
#           handle      ax
def plot1DSubgrid( name, grid, ax):
    length= grid.shape[0]
    x = np.arange(grid.shape[0])
    ax.plot(x,grid)
    # ax.set_title(str(grid.shape[0]), 'right')
    ax.set_xlim([0,length-1])
    ax.text(length/5, 0, " Max = " + str(np.round(np.max(grid))) + "\n Min = " +
        str(np.round(np.min(grid))) + "\n L=" + str(length))
    ax.locator_params(axis='y',nbins=16)
    ax.locator_params(axis='x',nbins=16)
    ax.grid()

    print name + "\t=" + str(np.max(grid))
