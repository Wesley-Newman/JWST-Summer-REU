import numpy as np
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
import stpsf
import matplotlib



def pre_generated_PSF_coord_fitting(fits_file, hdu, header, x_coords, y_coords, hdu_coarse_psf, hdu_fine_psf):
    """
    This function is used to fit the position of a star using a pre-generated point spread function. Given an initial guess for position of a star, this function with sistamatically check for the best fit for the position. 
    Parameters:
        fits_file : this is the fits file of the stars you are fitting. It is used to validate the answers by re-simulating a psf which pulls insturment data from the fits file

        hdu : this is the hdu that is being used

        header : header for the data used for x and y axis scale

        x_coords : this is the inital guess for the x coord of the star, this is automatically rounded to an interger but can be a float when calling function

        y_coords : this is the inital guess for the y coord of the star, this is automatically rounded to an interger but can be a float when calling function

        hdu_coarse_psf : this is the hdu for the coarse psf fitting. these psf's are pre-rendered to save compute time. View psf_generation.ipynb for info on the generation of these PSF's
    
        hdu_fine_psf : this is the hdu for the fine psf fitting. These psf's are pre-rendered to save compute time. View psf_generation.ipynb for info on the generation of these PSF's

    Outputs: 
        fine_best_x : this is the solution from the fine psf fitting solution. This is the best guess for the x coord of the star

        fine_best_y : this is the solution from the fine psf fitting solution. This is the best guess for the y coord of the star

    """
    #this defines the size of the cutout around the star that you are fitting and needs to match the size of the pre-rendered psfs
    boxsize = 25

    def psf_fitting(x_coords, y_coords, psf_library):
        #this rounding ensure that the same pixel is used evey time when a float is passed into Cutout2D
        x_coords = round(x_coords)
        y_coords = round(y_coords)
        psf_center = (x_coords, y_coords)

        #defines the region that is used for fitting
        obs_psf = Cutout2D(hdu['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
        obs_psf_err = Cutout2D(hdu['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()

        #background extraction and elimanation of nan and inf
        obs_psf -= np.nanmedian(obs_psf)          
        valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
        err_sq  = obs_psf_err[valid] ** 2     

        chi_squared = []

        x_shifts = []
        y_shifts = []

        #this loops through all of the psf that were generated and calculates the chi squared value of each shift
        for k in range(len(psf_library) - 1):
            x_shifts.append(psf_library[k+1].header['XSHIFT'])
            y_shifts.append(psf_library[k+1].header['ySHIFT'])

            sim_psf = psf_library[k+1].data
            sim_psf -= np.nanmedian(sim_psf)

            scalefactor = np.nansum(obs_psf) / np.nansum(sim_psf)


            diff = (obs_psf - (sim_psf * scalefactor))[valid]        

            chi_squared.append(np.nansum(diff**2/err_sq)/np.isfinite(diff).sum())

        #after all the chi squared values are computed we find the minimum
        min_value = np.argmin(chi_squared)
        min_chi = chi_squared[min_value]
        best_x = x_shifts[min_value]
        best_y = y_shifts[min_value]

        #defines the pixels per arcsecond along the x and y axis
        x_axis_scale = np.sqrt(header['CD1_1'] **2 + header['CD2_1'] **2 ) * 3600
        y_axis_scale = np.sqrt(header['CD1_2'] **2 + header['CD2_2'] **2 ) * 3600

        #this finds the new guess for the x and y value of the star in pixels
        x_pixel_shift = best_x / x_axis_scale
        y_pixel_shift = best_y / y_axis_scale
        new_x = x_coords + x_pixel_shift
        new_y = y_coords + y_pixel_shift

        #this code is for generating a contour plot to ensure that a global minimum is found
        #this reshapes the data into a 2d array
        chi_squared = np.array(chi_squared)
        n_grid = int(np.sqrt(len(psf_library) - 1))
        chi_squ_2d = np.array((chi_squared/ np.min(chi_squared))).reshape(n_grid, n_grid)
        x_grid_2d = np.array(x_shifts).reshape(n_grid, n_grid)
        y_grid_2d = np.array(y_shifts).reshape(n_grid, n_grid)

        #this plots the contour plot
        fig, ax = plt.subplots(figsize=(16, 9))
        levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13 ,14, 15, 16, 17,18, 19, 20]
        ct = ax.contour(x_grid_2d, y_grid_2d, chi_squ_2d, levels=levels, origin='lower')
        fig.colorbar(ct, ax=ax)
        ax.scatter(best_x, best_y, color='red', s=50,  zorder=1000, label=f'Δx={best_x}"  Δy={best_y}"')
        ax.set_xlabel('X offset (arcsec)')
        ax.set_ylabel('Y offset (arcsec)')
        ax.legend()
        plt.rc('font', size=10)
        plt.tight_layout()
        plt.show()
    

        return (new_x, new_y, min_chi, best_x, best_y, min_value)



    
    fine_best_x = []
    fine_best_y = []
    fine_x_shift = []
    fine_y_shift = []

    #when this function is called the prevous fitting code will be run for possibly 3 times to ensure a minimum is found
    #this loops through the prevous code and fits the position of the stars 
    for k in range(len(x_coords)):

        current_x = x_coords[k]
        current_y = y_coords[k]
        
        max_iter = 3

        for iteration in range(max_iter):

            coarse_x_fit, coarse_y_fit, _, _, _, _ = psf_fitting(current_x, current_y, hdu_coarse_psf)
     


            fine_x_fit, fine_y_fit, fine_min_chi, x_shift, y_shift, _ = psf_fitting(coarse_x_fit,coarse_y_fit,hdu_fine_psf)

            if fine_min_chi <= 50:

                fine_best_x.append(fine_x_fit)
                fine_best_y.append(fine_y_fit)
                fine_x_shift.append(x_shift)
                fine_y_shift.append(y_shift)

                break

            current_x = fine_x_fit
            current_y = fine_y_fit
        else:

            fine_best_x.append(fine_x_fit)
            fine_best_y.append(fine_y_fit)
            fine_x_shift.append(x_shift)
            fine_y_shift.append(y_shift)


    validation_chi = []
    #once a best position is found we can validate that solution through re-calculating the psf and subtracting it from the observed psf
    for k in range(len(fine_best_x)):
        rounded_x = round(fine_best_x[k])
        rounded_y = round(fine_best_y[k])

        psf_center = (rounded_x, rounded_y)

        inst = stpsf.setup_sim_to_match_file(fits_file, verbose= False)


        obs_psf = Cutout2D(hdu['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
        obs_psf_err = Cutout2D(hdu['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()

        obs_psf -= np.nanmedian(obs_psf)    # perform a simple background subtraction) 

        valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
        err_sq  = obs_psf_err[valid] ** 2    


        inst.options['source_offset_x'] = fine_x_shift[k]   
        inst.options['source_offset_y'] = fine_y_shift[k]

        sim_psf_hdu = inst.calc_psf(fov_pixels=boxsize)
        sim_psf = sim_psf_hdu['DET_DIST'].data
        sim_psf -= np.nanmedian(sim_psf)

        scalefactor = np.nansum(obs_psf) / np.nansum(sim_psf)

        diff = (obs_psf - (sim_psf * scalefactor))[valid]
        chisqr = np.nansum(diff**2/err_sq)/np.isfinite(diff).sum()


        diff_display = obs_psf - (sim_psf * scalefactor)

        fig, axes = plt.subplots(figsize=(13,3), ncols=3)

        vmax = np.nanmax(obs_psf)
        cmap = matplotlib.cm.gist_heat
        cmap.set_bad(cmap(0))
        axes[0].imshow(obs_psf, norm = matplotlib.colors.LogNorm(vmax/1e4, vmax), cmap=cmap, origin='lower')
        axes[0].set_title("Observed PSF from science data")

        stpsf.display_psf(sim_psf_hdu, ext='DET_DIST', vmax=0.1, vmin=1e-5, ax=axes[1], )

        axes[2].imshow(diff_display, norm = matplotlib.colors.LogNorm(vmax/1e5, vmax), cmap=cmap, origin='lower')


        validation_chi.append(chisqr)

        axes[2].set_title("Difference")
        axes[2].text(3, 3, f"$\\chi^2$ = {chisqr}", color='white')
        for ax in [axes[0], axes[2]]:
            ax.set_xlabel("Pixels")


    return(fine_best_x, fine_best_y)