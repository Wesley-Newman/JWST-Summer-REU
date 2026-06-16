import numpy as np
from astropy.io import fits
import stpsf
from astropy.nddata import Cutout2D
from tqdm import tqdm
import astropy.wcs.utils as utils
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib

def pre_generated_PSF_coord_fitting(fits_file, hdu, header, x_coords, y_coords, hdu_coarse_psf, hdu_fine_psf):
    def psf_fitting(x_coords, y_coords, psf_library):
        x_coords = np.round(x_coords)
        y_coords = np.round(y_coords)

        psf_center = (x_coords, y_coords)
        boxsize = 50

        obs_psf = Cutout2D(hdu['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
        obs_psf_err = Cutout2D(hdu['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()

        obs_psf -= np.nanmedian(obs_psf)          

        valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
        # n_valid = int(valid.sum())
        err_sq  = obs_psf_err[valid] ** 2     

        chi_squared = []

        x_shifts = []
        y_shifts = []

        for k in range(len(psf_library) - 1):
            x_shifts.append(psf_library[k+1].header['XSHIFT'])
            y_shifts.append(psf_library[k+1].header['ySHIFT'])

            sim_psf = psf_library[k+1].data

            scalefactor = np.nansum(obs_psf) / np.nansum(sim_psf)


            diff = (obs_psf - (sim_psf * scalefactor))[valid]        

            chi_squared.append(np.nansum(diff**2/err_sq)/np.isfinite(diff).sum())

        min_value = np.argmin(chi_squared)
        min_chi = chi_squared[min_value]
        best_x = x_shifts[min_value]
        best_y = y_shifts[min_value]



        x_axis_scale = np.sqrt(header['CD1_1'] **2 + header['CD2_1'] **2 ) * 3600
        y_axis_scale = np.sqrt(header['CD1_2'] **2 + header['CD2_2'] **2 ) * 3600


        x_pixel_shift = best_x / x_axis_scale
        y_pixel_shift = best_y / y_axis_scale

        new_x = x_coords + x_pixel_shift
        new_y = y_coords + y_pixel_shift


        # print (new_x, new_y)


        chi_squared = np.array(chi_squared)

        n_grid = int(np.sqrt(len(psf_library) - 1))

        chi_squ_2d = np.array((chi_squared/ np.min(chi_squared))).reshape(n_grid, n_grid)


        x_grid_2d = np.array(x_shifts).reshape(n_grid, n_grid)
        y_grid_2d = np.array(y_shifts).reshape(n_grid, n_grid)

        fig, ax = plt.subplots(figsize=(16, 9))

        levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13 ,14, 15, 16, 17,18, 19, 20]

        ct = ax.contour(x_grid_2d, y_grid_2d, chi_squ_2d, levels=levels, origin='lower')

        fig.colorbar(ct, ax=ax)
        ax.scatter(best_x,   best_y,   color='red',    s=50,  zorder=1000,
                   label=f'Δx={best_x}"  Δy={best_y}"')

        ax.set_xlabel('X offset (arcsec)')
        ax.set_ylabel('Y offset (arcsec)')
        ax.legend()
        plt.rc('font', size=10)
        plt.tight_layout()
        plt.show()




        # k_best = np.argmin(chi_squared)
        # print(f"Best k: {k_best}")
        # print(f"Best shift X: {psf_library[k_best+1].header['XSHIFT']}")
        # print(f"Best shift Y: {psf_library[k_best+1].header['ySHIFT']}")
        # print(f"fine_min_chi: {chi_squared[k_best]:.6f}")

    

        return (new_x, new_y, min_chi, best_x, best_y, min_value)


    best_x = []
    best_y = []
    min_chi = []

    
    fine_best_x = []
    fine_best_y = []
    fine_min_chi = []
    fine_x_shift = []
    fine_y_shift = []
    index = []

    for k in range(len(x_coords)):

        current_x = x_coords[k]
        current_y = y_coords[k]
        
        max_iter = 10

        for iteration in range(max_iter):

            bestll_x, bestll_y, min_chill, _, _, _ = psf_fitting(current_x, current_y, hdu_coarse_psf)
     


            bestl_x, bestl_y, min_chil, x_shiftl, y_shiftl, lindex = psf_fitting(
                bestll_x,
                bestll_y,
                hdu_fine_psf
            )

            # If fit is acceptable, save final result and exit loop
            if min_chil <= 50:
                best_x.append(bestll_x)
                best_y.append(bestll_y)
                min_chi.append(min_chill)
                fine_best_x.append(bestl_x)
                fine_best_y.append(bestl_y)
                fine_min_chi.append(min_chil)
                fine_x_shift.append(x_shiftl)
                fine_y_shift.append(y_shiftl)
                index.append(lindex)
                break

            # Otherwise use this result as the next starting point
            current_x = bestl_x
            current_y = bestl_y
        else:
            best_x.append(bestll_x)
            best_y.append(bestll_y)
            min_chi.append(min_chill)
            fine_best_x.append(bestl_x)
            fine_best_y.append(bestl_y)
            fine_min_chi.append(min_chil)
            fine_x_shift.append(x_shiftl)
            fine_y_shift.append(y_shiftl)
            index.append(lindex)        


    # validation_chi = []
    
    # for k in range(len(best_x)):
    #     rounded_x = np.round(best_x[k])
    #     rounded_y = np.round(best_y[k])

    #     # x_offset_pix = fine_best_x[k] - rounded_x
    #     # y_offset_pix = fine_best_y[k] - rounded_y

    #     # x_axis_scale = np.sqrt(header['CD1_1'] **2 + header['CD2_1'] **2 ) * 3600
    #     # y_axis_scale = np.sqrt(header['CD1_2'] **2 + header['CD2_2'] **2 ) * 3600

    #     # x_offset_arc = x_offset_pix * x_axis_scale
    #     # y_offset_arc = y_offset_pix * y_axis_scale.copy()

    #     psf_center = (rounded_x, rounded_y)

    #     boxsize = 50
    #     inst = stpsf.setup_sim_to_match_file(fits_file, verbose= False)

    #     # Load that science data. Cut out the surface brightness and uncertainty for the desired location
    #     # obs_im = fits.open(fits_file)

    #     obs_psf = Cutout2D(hdu['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
    #     obs_psf_err = Cutout2D(hdu['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()

    #     obs_psf -= np.nanmedian(obs_psf)    # perform a simple background subtraction) 

    #     valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
    #     n_valid = int(valid.sum())
    #     err_sq  = obs_psf_err[valid] ** 2    

    #     # print (fine_x_shift, fine_y_shift)

    #     inst.options['source_offset_x'] = fine_x_shift[k]   # in units of arcseconds. So, this is about 1/5 of a NRC LW pixel
    #     inst.options['source_offset_y'] = fine_y_shift[k]

    #     sim_psf = inst.calc_psf(fov_pixels=boxsize)

    #     scalefactor = np.nansum(obs_psf) / np.nansum(sim_psf['DET_DIST'].data)

    #     diff = (obs_psf - (sim_psf['DET_DIST'].data * scalefactor))[valid]
    #     chisqr = np.nansum(diff**2/err_sq)/np.isfinite(diff).sum()

    #     diff_display = obs_psf - (sim_psf['DET_DIST'].data * scalefactor)

    #     fig, axes = plt.subplots(figsize=(13,3), ncols=3)

    #     vmax = np.nanmax(obs_psf)
    #     cmap = matplotlib.cm.gist_heat
    #     cmap.set_bad(cmap(0))
    #     axes[0].imshow(obs_psf, norm = matplotlib.colors.LogNorm(vmax/1e4, vmax), cmap=cmap, origin='lower')
    #     axes[0].set_title("Observed PSF from science data")

    #     stpsf.display_psf(sim_psf, ext='DET_DIST', vmax=0.1, vmin=1e-5, ax=axes[1], )

    #     axes[2].imshow(diff_display, norm = matplotlib.colors.LogNorm(vmax/1e5, vmax), cmap=cmap, origin='lower')


    #     validation_chi.append(chisqr)

    #     axes[2].set_title("Difference")
    #     axes[2].text(3, 3, f"$\\chi^2$ = {chisqr}", color='white')
    #     for ax in [axes[0], axes[2]]:
    #         ax.set_xlabel("Pixels")

    # if np.all(np.array(validation_chi) == np.array(fine_min_chi)):
    #     return
    # else:
    #     print('AHHHHHHHHHHHHHHHHHH BAD BAD BAD THIS IS HORRIBLE AHHHHHHHHHH')
    return(fine_best_x, fine_best_y)