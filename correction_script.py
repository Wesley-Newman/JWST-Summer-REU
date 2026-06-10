import numpy as np
import astropy.units as u
import astropy.wcs.utils as utils
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
# from astroquery.gaia import Gaia
# Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib
import stpsf
from scipy.ndimage import shift as ndimage_shift
import warnings
from astropy.wcs import FITSFixedWarning
warnings.simplefilter('ignore', FITSFixedWarning)
from astropy.wcs import WCS
from astropy.time import Time

# def gaia_find_stars(ra_center, dec_center, size_deg, image, header):

#     coords = SkyCoord(ra = ra_center, dec = dec_center, unit=(u.degree, u.degree), frame='icrs')

#     gaia_star_field = Gaia.query_object_async(coordinate=coords, width= u.Quantity(size_deg, u.deg), height= u.Quantity(size_deg, u.deg))

#     skycoord_object = SkyCoord(ra = gaia_star_field['ra'], dec = gaia_star_field['dec'], unit=(u.degree, u.degree), frame='icrs')

#     image_gaia_stars = utils.skycoord_to_pixel(skycoord_object, WCS(header))

#     x = image_gaia_stars[0].astype(int)
#     y = image_gaia_stars[1].astype(int)

#     mask = (
#         (x >= 0) & (x < image.shape[1]) &
#         (y >= 0) & (y < image.shape[0])
#     )

#     mask[mask] &= ~np.isnan(image[y[mask], x[mask]])



#     gaia_star_field = gaia_star_field[mask]
#     x = x[mask]
#     y = y[mask]

#     fig, ax = plt.subplots(1, 1, figsize=(9, 9))
#     ax.imshow(image, origin='lower', vmin=np.nanmin(image), vmax=np.nanmax(image)*.01, cmap='viridis', interpolation='none')

#     ax.scatter(x, y, color = 'red')

#     plt.xticks([])
#     plt.yticks([])
#     ax.spines['top'].set_visible(False)
#     ax.spines['bottom'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.spines['left'].set_visible(False)
#     plt.rc('font', size=10)

#     plt.show()

#     return x, y, gaia_star_field


def find_stars(x_coords, y_coords, image):


    fig, ax = plt.subplots(1, 1, figsize=(9, 9))
    ax.imshow(image, origin='lower', vmin=np.nanmin(image), vmax=np.nanmax(image)*.01, cmap='viridis', interpolation='none')

    ax.scatter(x_coords, y_coords, color = 'red')

    plt.xticks([])
    plt.yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.rc('font', size=10)

    plt.show()


# def apply_proper_motion_correction(header, gaia_star):

#     GAIA_DR3_EPOCH = 2016.0                                    
#     obs_epoch = Time(header['MJD-AVG'], format='mjd').jyear   
#     dt = obs_epoch - GAIA_DR3_EPOCH                            

#     ra_gaia  = gaia_star['ra']   
#     dec_gaia = gaia_star['dec']   
#     pmra     = gaia_star['pmra'] 
#     pmdec    = gaia_star['pmdec']

    
#     dec_rad  = np.deg2rad(dec_gaia)
#     ra_corr  = ra_gaia  + (pmra  * dt) / (3600000 * np.cos(dec_rad))  
#     dec_corr = dec_gaia + (pmdec * dt) /  3600000                      

#     wcs = WCS(header)
#     x_pm, y_pm = wcs.all_world2pix(ra_corr, dec_corr, 0) 

#     return float(x_pm), float(y_pm)


# def gaia_star_check(image, header, gaia_star):


#     for index in range(len(gaia_star)):

#         fig, ax = plt.subplots(1, 1,figsize=(9, 9))

#         ax.imshow(image, origin='lower', vmin=0, vmax=100, cmap='viridis', interpolation='none')
        
#         x_pm, y_pm = apply_proper_motion_correction(header, gaia_star[index])

#         ax.set_xlim([x_pm - 25, x_pm + 25])
#         ax.set_ylim([y_pm - 25, y_pm + 25])


#         ax.scatter(x_pm, + y_pm, color='red')
#         ax.spines['top'].set_visible(False)
#         ax.spines['bottom'].set_visible(False)
#         ax.spines['right'].set_visible(False)
#         ax.spines['left'].set_visible(False)

#         plt.rc('font', size=10)

#         plt.show()


def star_check(image, x_coords, y_coords):

    x_round = np.round(x_coords)
    y_round = np.round(y_coords)

    size = 5
    
    for k in range(len(x_coords)):

        fig, ax = plt.subplots(1, 1,figsize=(9, 9))

        ax.imshow(image, origin='lower', vmin=0, vmax=2000, cmap='viridis', interpolation='none')
        plt.grid(visible=None, which='major', axis='both')

        
        # x_pm, y_pm = apply_proper_motion_correction(header, gaia_star[index])

        ax.set_xlim([x_round[k] - size, x_round[k] + size])
        ax.set_ylim([y_round[k] - size, y_round[k] + size])


        ax.scatter(x_coords[k], y_coords[k], color='red', marker='+', s=100, linewidths=2)
        # ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # ax.spines['left'].set_visible(False)

        plt.rc('font', size=10)
        plt.show()


# def gaia_psf_fitting(fits_file, header, shift_amount, test_size, gaia_star_field):
    

 
    # inst = stpsf.setup_sim_to_match_file(fits_file, verbose=False)

    # if gaia_star_field is not None:
    #     x_coord, y_coord = apply_proper_motion_correction(header, gaia_star_field)

    # psf_center = (x_coord, y_coord)
    # boxsize = 50
 
    # with fits.open(fits_file) as hdul:
    #     obs_psf     = Cutout2D(hdul['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
    #     obs_psf_err = Cutout2D(hdul['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
 
    # obs_psf    -= np.nanmedian(obs_psf)          
 
    # valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
    # n_valid = int(valid.sum())
    # err_sq  = obs_psf_err[valid] ** 2           
 
    # DET_MAX = 2047
    # inst.detector_position = (
    #     int(np.clip(x_coord, 0, DET_MAX)),
    #     int(np.clip(y_coord, 0, DET_MAX))
    # )

    # inst.options['source_offset_x'] = 0.0
    # inst.options['source_offset_y'] = 0.0
 
    # base_hdul   = inst.calc_psf(fov_pixels=boxsize)
    # base_psf    = base_hdul["DET_DIST"].data
    # pixel_scale = header['CDELT1']*3600    # arcsec / pixel


    # def chi_sq(dx_arcsec: float, dy_arcsec: float) -> float:
        
    #     d_col =  dx_arcsec / pixel_scale  
    #     d_row =  dy_arcsec / pixel_scale  
    #     shifted = ndimage_shift(base_psf, [d_row, d_col],
    #                             order=3, mode='constant', cval=0.0)


    #     scalefactor = np.nansum(obs_psf) / np.nansum(shifted)
        
    #     diff = (obs_psf - shifted * scalefactor)[valid]
    #     return float(np.nansum(diff**2 / err_sq) / n_valid)

    
    # shift_values = []

    # for p in range(-test_size, test_size + 1, +1):
    #     for j in range(test_size, -test_size - 1, -1):
    #         shift_values.append([j*shift_amount, p*shift_amount])
 
    # chi_squared = [
    #     chi_sq(sv[1], sv[0])                         
    #     for sv in tqdm(shift_values)
    # ]

    # min_value = np.argwhere(chi_squared == np.min(chi_squared))


    # index = min_value[0, 0]

    # chi_squared = np.array(chi_squared)

    # chi_squ_2d = np.array((chi_squared/ np.min(chi_squared))).reshape((2 * test_size + 1, 2 * test_size + 1))

    # x_grid = []
    # y_grid = []

    # for k, ds in enumerate(shift_values):
    #     x_grid.append(ds[1])
    #     y_grid.append(ds[0])

    # x_grid = np.array(x_grid)
    # y_grid = np.array(y_grid)

    # x_grid_2d = np.array(x_grid).reshape((2 * test_size + 1, 2 * test_size + 1))
    # y_grid_2d = np.array(y_grid).reshape((2 * test_size + 1, 2 * test_size + 1))

    # vx_left = np.min(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    # vx_right = np.max(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    # vy_bottom = np.min(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    # vy_top = np.max(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])

    # h_error = np.abs(vx_left - vx_right )
    # w_error = np.abs(vy_bottom - vy_top)

    # error_area = h_error * w_error
 
    # fig, ax = plt.subplots(figsize=(16, 9))


    # levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13 ,14, 15, 16, 17,18, 19, 20]

    # ct = ax.contour(x_grid_2d, y_grid_2d, chi_squ_2d, levels=levels, origin='lower')

    # fig.colorbar(ct, ax=ax)
    # ax.scatter(shift_values[index][1],   shift_values[index][0],   color='red',    s=50,  zorder=1000,
    #            label=f'Δx={shift_values[index][1]:.4f}"  Δy={shift_values[index][0]:.4f}"')
    
    # ax.set_xlabel('X offset (arcsec)')
    # ax.set_ylabel('Y offset (arcsec)')
    # ax.legend()
    # plt.rc('font', size=10)
    # plt.tight_layout()
    # plt.show()
    
    # chi_squared_min = np.min(chi_squared)
    # return chi_squared_min, shift_values[index][1], shift_values[index][0], error_area, h_error, w_error

def psf_fitting(fits_file, header, shift_amount, test_size, x_coord, y_coord):
    

    # while True
    #     if not np.any(outliers):
    #         break

    # mask[np.where(mask)[0][outliers]] = False


    inst = stpsf.setup_sim_to_match_file(fits_file, verbose=False)

    x_coord = np.round(x_coord)
    y_coord = np.round(y_coord)

    psf_center = (x_coord, y_coord)

    center_original = utils.pixel_to_skycoord(x_coord, y_coord, WCS(header), origin=0, mode='all', cls=None)

    boxsize = 50
 
    with fits.open(fits_file) as hdul:
        sci_cutout  = Cutout2D(hdul['SCI'].data, position=psf_center, size=boxsize, wcs=WCS(header), mode='partial', fill_value=np.nan)
        obs_psf     = sci_cutout.data.copy()
        obs_psf_wcs = sci_cutout.wcs          
        obs_psf_err = Cutout2D(hdul['ERR'].data, position=psf_center, size=boxsize, wcs=WCS(header), mode='partial', fill_value=np.nan).data.copy()

    

    obs_psf -= np.nanmedian(obs_psf)          

    fig, ax = plt.subplots(1, 1,figsize=(9, 9))

    ax.imshow(obs_psf, origin='lower', vmin=0, vmax=2000, cmap='viridis', interpolation='none')
    plt.grid(visible=None, which='major', axis='both')

 
    valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
    n_valid = int(valid.sum())
    err_sq  = obs_psf_err[valid] ** 2           
 

    inst.options['source_offset_x'] = 0.0
    inst.options['source_offset_y'] = 0.0
 
    base_hdul   = inst.calc_psf(fov_pixels=boxsize)
    base_psf    = base_hdul["DET_DIST"].data

    def chi_sq(dx: float, dy: float) -> float:
        
        d_col =  dx  
        d_row =  dy  
        shifted = ndimage_shift(base_psf, [d_row, d_col],
                                order=3, mode='constant', cval=0.0)


        scalefactor = np.nansum(obs_psf) / np.nansum(shifted)
        
        diff = (obs_psf - shifted * scalefactor)[valid]
        return float(np.nansum(diff**2 / err_sq) / n_valid)

    
    shift_values = []

    for p in range(-test_size, test_size + 1, +1):
        for j in range(test_size, -test_size - 1, -1):
            shift_values.append([j*shift_amount, p*shift_amount])
 
    chi_squared = [
        chi_sq(sv[1], sv[0])                         
        for sv in tqdm(shift_values)
    ]

    min_value = np.argwhere(chi_squared == np.min(chi_squared))


    index = min_value[0, 0]

    chi_squared = np.array(chi_squared)

    chi_squ_2d = np.array((chi_squared/ np.min(chi_squared))).reshape((2 * test_size + 1, 2 * test_size + 1))

    x_grid = []
    y_grid = []

    for k, ds in enumerate(shift_values):
        x_grid.append(ds[1])
        y_grid.append(ds[0])

    x_grid = np.array(x_grid)
    y_grid = np.array(y_grid)

    x_grid_2d = np.array(x_grid).reshape((2 * test_size + 1, 2 * test_size + 1))
    y_grid_2d = np.array(y_grid).reshape((2 * test_size + 1, 2 * test_size + 1))

    vx_left = np.min(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    vx_right = np.max(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    vy_bottom = np.min(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
    vy_top = np.max(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])

    h_error = np.abs(vx_left - vx_right )
    w_error = np.abs(vy_bottom - vy_top)

    error_area = h_error * w_error
 
    fig, ax = plt.subplots(figsize=(16, 9))


    levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13 ,14, 15, 16, 17,18, 19, 20]

    ct = ax.contour(x_grid_2d, y_grid_2d, chi_squ_2d, levels=levels, origin='lower')

    fig.colorbar(ct, ax=ax)
    ax.scatter(shift_values[index][1],   shift_values[index][0],   color='red',    s=50,  zorder=1000,
               label=f'Δx={shift_values[index][1]:.2f}"  Δy={shift_values[index][0]:.2f}"')
    
    ax.set_xlabel('X offset (arcsec)')
    ax.set_ylabel('Y offset (arcsec)')
    ax.legend()
    plt.rc('font', size=10)
    plt.tight_layout()
    plt.show()
    
    chi_squared_min = np.min(chi_squared)


    origin = utils.pixel_to_skycoord(0.0, 0.0, obs_psf_wcs, origin=0, mode='all', cls=None)
    
    found_shift = utils.pixel_to_skycoord(0.0 + shift_values[index][1], 0.0 + shift_values[index][0], obs_psf_wcs, origin=0, mode='all', cls=None)

    dif_ra = (found_shift.ra.deg - origin.ra.deg)
    dif_dec = (found_shift.dec.deg - origin.dec.deg)

    new_ra = center_original.ra.deg + dif_ra
    new_dec = center_original.dec.deg + dif_dec

    new_skycoord = SkyCoord(ra = new_ra, dec = new_dec, unit = 'deg')
    
    new_x_coord, new_y_coord = utils.skycoord_to_pixel(new_skycoord, WCS(header))

    # x_pixel_shift = shift_values[index][1] #- initial[1]
    # y_pixel_shift = shift_values[index][0] #- initial[0]

    # new_x_coord = x_pixel_shift + x_coord
    # new_y_coord = y_pixel_shift + y_coord

    # print(f"Original: ({x_coord:.0f}, {y_coord:.0f}) "
    # f"Shift: ({x_pixel_shift:.2f}, {y_pixel_shift:.2f})")

    return chi_squared_min, new_x_coord, new_y_coord, error_area, h_error, w_error




# def better_psf_fitting(fits_file, header, shift_amount, test_size, x_coord, y_coord):
    

 
#     inst = stpsf.setup_sim_to_match_file(fits_file, verbose=False)

#     psf_center = (x_coord, y_coord)
#     boxsize = 50
 
#     with fits.open(fits_file) as hdul:
#         cutout_sci = Cutout2D(hdul['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan)
#         cutout_err = Cutout2D(hdul['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan)   
#         obs_psf     = Cutout2D(hdul['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()
#         obs_psf_err = Cutout2D(hdul['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data.copy()

#         actual_center_x = cutout_sci.center_original[0]
#         actual_center_y = cutout_sci.center_original[1]
 
#     obs_psf    -= np.nanmedian(obs_psf)            

#     inst.options['source_offset_x'] = 0.0
#     inst.options['source_offset_y'] = 0.0
 
#     base_hdul   = inst.calc_psf(fov_pixels=boxsize)
#     base_psf    = base_hdul["DET_DIST"].data

#     pixel_scale = header['CDELT1'] * 3600


#     # obs_psf = obs_psf[2:-2, 2:-2]
#     # obs_psf_err = obs_psf_err[2:-2, 2:-2]

#     valid   = np.isfinite(obs_psf) & np.isfinite(obs_psf_err) & (obs_psf_err > 0)
#     n_valid = int(valid.sum())
#     err_sq  = obs_psf_err[valid] ** 2         

#     def chi_sq(dx_arcsec, dy_arcsec):
        
#         d_x =  dx_arcsec / pixel_scale  
#         d_y =  dy_arcsec / pixel_scale  

#         shifted = ndimage_shift(base_psf, [d_y, d_x],
#                                 order=3, mode='constant', cval=0.0)

#         # shifted = shifted[2:-2, 2:-2]

#         scalefactor = np.nansum(obs_psf[valid]) / np.nansum(shifted[valid])
        
#         diff = (obs_psf - shifted * scalefactor)[valid]
#         return float(np.nansum(diff**2 / err_sq) / n_valid)

    
#     shift_values = []

#     for p in range(-test_size, test_size + 1, +1):
#         for j in range(test_size, -test_size - 1, -1):
#             shift_values.append([j*shift_amount, p*shift_amount])
 
#     chi_squared = [
#         chi_sq(sv[1], sv[0])                         
#         for sv in tqdm(shift_values)
#     ]

#     min_value = np.argwhere(chi_squared == np.min(chi_squared))

#     index = min_value[0, 0]

#     chi_squared = np.array(chi_squared)

#     chi_squ_2d = np.array((chi_squared/ np.min(chi_squared))).reshape((2 * test_size + 1, 2 * test_size + 1))

#     x_grid = []
#     y_grid = []

#     for k, ds in enumerate(shift_values):
#         x_grid.append(ds[1])
#         y_grid.append(ds[0])

#     x_grid = np.array(x_grid)
#     y_grid = np.array(y_grid)

#     x_grid_2d = np.array(x_grid).reshape((2 * test_size + 1, 2 * test_size + 1))
#     y_grid_2d = np.array(y_grid).reshape((2 * test_size + 1, 2 * test_size + 1))

#     vx_left = np.min(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
#     vx_right = np.max(x_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
#     vy_bottom = np.min(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])
#     vy_top = np.max(y_grid[np.argwhere((chi_squared/(np.min(chi_squared)) <= 2))])

#     h_error = np.abs(vx_left - vx_right )
#     w_error = np.abs(vy_bottom - vy_top)

#     error_area = h_error * w_error
 
#     fig, ax = plt.subplots(figsize=(16, 9))


#     levels = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ,13 ,14, 15, 16, 17,18, 19, 20]

#     ct = ax.contour(x_grid_2d, y_grid_2d, chi_squ_2d, levels=levels, origin='lower')

#     fig.colorbar(ct, ax=ax)
#     ax.scatter(shift_values[index][1],   shift_values[index][0],   color='red',    s=50,  zorder=1000,
#                label=f'Δx={shift_values[index][1]:.4f}"  Δy={shift_values[index][0]:.4f}"')
    
#     ax.set_xlabel('X offset (arcsec)')
#     ax.set_ylabel('Y offset (arcsec)')
#     ax.legend()
#     plt.rc('font', size=10)
#     plt.tight_layout()
#     plt.show()
    
#     chi_squared_min = np.min(chi_squared)

#     x_pix_shift = shift_values[index][1] / pixel_scale
#     y_pix_shift = shift_values[index][0] / pixel_scale

#     x_true = actual_center_x + x_pix_shift
#     y_true = actual_center_y + y_pix_shift
#     print(
#     f"Original: ({x_coord:.3f}, {y_coord:.3f}) "
#     f"Shift: ({x_pix_shift:.3f}, {y_pix_shift:.3f})")
#     return chi_squared_min, (x_true), (y_true), error_area, h_error, w_error

#     # return chi_squared_min, (x_coord + x_pix_shift), (y_coord + y_pix_shift), error_area, h_error, w_error


def sigma_clip(x, *arrays, threshold=3):

    x = np.asarray(x, dtype=float)
    mask = np.ones(len(x), dtype=bool)

    while True:
        current = x[mask]
        z_scores = np.abs((current - np.mean(current)) / np.std(current))

        outliers = z_scores > threshold

        if not np.any(outliers):
            break

        mask[np.where(mask)[0][outliers]] = False

    return x[mask], *[arr[mask] for arr in arrays]
 


# def gaia_plot_data_sim_comparison(masked_x_coords, masked_y_coords, file, header, gaia_star_field, masked_best_x, masked_best_y):

#     if gaia_star_field is not None:
#         x_coord, y_coord = apply_proper_motion_correction(header, gaia_star_field)


#     psf_center = (masked_x_coords, masked_y_coords)
#     boxsize = 50
#     inst = stpsf.setup_sim_to_match_file(file, verbose= False)

#     # Load that science data. Cut out the surface brightness and uncertainty for the desired location
#     obs_im = fits.open(file)

#     obs_psf = Cutout2D(obs_im['SCI'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data
#     obs_psf_err = Cutout2D(obs_im['ERR'].data, position=psf_center, size=boxsize, mode='partial', fill_value=np.nan).data

#     obs_psf -= np.nanmedian(obs_psf)    # perform a simple background subtraction) 

#     inst.options['source_offset_x'] = masked_best_x   # in units of arcseconds. So, this is about 1/5 of a NRC LW pixel
#     inst.options['source_offset_y'] = masked_best_y
#     sim_psf_offset = inst.calc_psf(fov_pixels=boxsize)

#     fig, axes = plt.subplots(figsize=(13,3), ncols=3)

#     vmax = np.nanmax(obs_psf)
#     cmap = matplotlib.cm.gist_heat
#     cmap.set_bad(cmap(0))
#     axes[0].imshow(obs_psf, norm = matplotlib.colors.LogNorm(vmax/1e4, vmax), cmap=cmap, origin='lower')
#     axes[0].set_title("Observed PSF from science data")

#     stpsf.display_psf(sim_psf_offset, ext='DET_DIST', vmax=0.1, vmin=1e-5, ax=axes[1], )

#     scalefactor = np.nansum(obs_psf) / np.sum(sim_psf_offset["DET_DIST"].data)
#     difference = obs_psf - sim_psf_offset["DET_DIST"].data *scalefactor
#     axes[2].imshow(difference, norm = matplotlib.colors.LogNorm(vmax/1e5, vmax), cmap=cmap, origin='lower')

#     chisqr = np.nansum(difference**2/obs_psf_err**2)/np.isfinite(difference).sum()
#     axes[2].set_title("Difference")
#     axes[2].text(3, 3, f"$\\chi^2$ = {chisqr:.2f}", color='white')
#     for ax in [axes[0], axes[2]]:
#         ax.set_xlabel("Pixels")

def plot_data_sim_comparison(obs_psf, obs_psf_err, sim_psf_offset):
    fig, axes = plt.subplots(figsize=(13,3), ncols=3)

    vmax = np.nanmax(obs_psf)
    cmap = matplotlib.cm.gist_heat
    cmap.set_bad(cmap(0))
    axes[0].imshow(obs_psf, norm = matplotlib.colors.LogNorm(vmax/1e4, vmax), cmap=cmap, origin='lower')
    axes[0].set_title("Observed PSF from science data")

    stpsf.display_psf(sim_psf_offset, ext='DET_DIST', vmax=0.1, vmin=1e-5, ax=axes[1], )

    scalefactor = np.nansum(obs_psf) / np.sum(sim_psf_offset["DET_DIST"].data)
    difference = obs_psf - sim_psf_offset["DET_DIST"].data *scalefactor
    axes[2].imshow(difference, norm = matplotlib.colors.LogNorm(vmax/1e5, vmax), cmap=cmap, origin='lower')

    chisqr = np.nansum(difference**2/obs_psf_err**2)/np.isfinite(difference).sum()
    axes[2].set_title("Difference")
    axes[2].text(3, 3, f"$\\chi^2$ = {chisqr:.2f}", color='white')
    for ax in [axes[0], axes[2]]:
        ax.set_xlabel("Pixels")

