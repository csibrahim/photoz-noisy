pro plot_galaxy, catalog, id, _extra=extra

    if catalog eq 'cosmos' then catalog = 'COSMOS_data_PH_27_11_18'
    if catalog eq 'uds'    then catalog =  'XMM_data_PH_27_11_18'

    result_dir = '../results/eazy_'+catalog+'/'

    readcol, result_dir+'seds/eazy_'+strn(id)+'.temp_sed', lam, sed, lam0, sed0
    readcol, result_dir+'seds/eazy_'+strn(id)+'.obs_sed', blam, bflx, berr, bfull, bmod, bmod0

    plot, lam*1d-4, ujy2cgs(lam, sed)*1d19, /xlog, xtit='wavelength [um]', $
        yr=[-2*median(ujy2cgs(blam, bfull))*1d19, max(ujy2cgs(blam,bflx))*1.2*1d19]

    errplot, blam*1d-4, ujy2cgs(blam, bflx-bfull)*1d19, ujy2cgs(blam, bflx+bfull)*1d19, col='ff'x
    oplot, blam*1d-4, ujy2cgs(blam, bflx)*1d19, psym=symcat(16), col='ff'x

    oplot, blam*1d-4, ujy2cgs(blam, bmod)*1d19, psym=symcat(16), col='ff0000'x
end
