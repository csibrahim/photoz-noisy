catalogs=['COSMOS_data_PH_27_11_18', 'XMM_data_PH_27_11_18']

zbins=[0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.500]
nbin = n_elements(zbins)-1

for c=0, n_elements(catalogs)-1 do begin
    catalog = catalogs[c]
    result_dir = '../results/eazy_'+catalog+'/'

    zphot_file = result_dir+'eazy.zout'

    if ~file_test(zphot_file) then begin
        print, 'could not find ', zphot_file
        continue
    endif

    zphot = readtable(zphot_file, $
        types=['L',replicate('F',14),'I',replicate('F',4)], $
        names=['id','z_spec','z_a','z_m1','chi_a','z_p','chi_p','z_m2','odds','l68','u68','l95','u95','l99','u99','nfilt','q_z','z_peak','peak_prob','z_mc'])

    znames = ['z_a', 'z_m2', 'z_peak']
    zshow = ['z_chi2', 'z_mean', 'z_mode']

    ; z_a    : best chi2
    ; z_m1   : p(z) weighted (no prior)
    ; z_p    : best chi2 + prior
    ; z_m2   : p(z) weighted (with prior)
    ; z_peak : z of peak p(z)

    out_stack_dir = result_dir+'stacked_pz/'
    file_mkdir, out_stack_dir

    for iz=0, n_elements(znames)-1 do begin
        icz = (where(tag_names(zphot) eq strupcase(znames[iz]), cnt))[0]
        if cnt eq 0 then message, 'could not get column '+znames[iz]

        zprint = strreplace(zshow[iz], '_', '!d')
        z = zphot.(icz)

        zmax = 5.0
        zlim = 9.0
        dz = 0.01

        while psplot(result_dir+'zspec_'+zshow[iz]+'.eps', [8.5,8]) do begin
            !x.margin = [7,3]

            densplot, z, zphot.z_spec, xr=[0,zmax], yr=[0,zmax], numbins=zmax*50, $
                xtit=zprint, ytit='z!dspec', /log, $
                coltable=['88ff'x, !p.color], xdata=xdata, nan='ffffff'x, $
                /zero_as_nan, weight=replicate(1.0, n_elements(z))

            oplotline, 1, 0, thick=2, col='ff'x
        endwhile

        np = 5
        while psplot(result_dir+'metrics_'+zshow[iz]+'.eps', [8.5,4*np]) do begin
            !y.margin = [0,0]
            !y.omargin = [4,1]

            pos = mplotpos(multi=[1,np])

            zx = fltarr(2*nbin)
            zx[2*indgen(nbin)+0] = zbins[0:nbin-1]
            zx[2*indgen(nbin)+1] = zbins[1:nbin]

            for p=0, np-1 do begin
                if p eq np-1 then xtit = zprint else xtit = ''
                if p eq np-1 then xtickformat = '' else xtickformat = 'no_tick'

                if p eq 0 then ytit = 'bias'
                if p eq 1 then ytit = 'scatter'
                if p eq 2 then ytit = '% outlier'
                if p eq 3 then ytit = '% of |'+cggreek('Delta')+'z| < '+cggreek('sigma')+'!dz!n(68%)'
                if p eq 4 then ytit = '% of |'+cggreek('Delta')+'z| < '+cggreek('sigma')+'!dz!n(90%)'

                my = fltarr(2*nbin)
                for ib=0, nbin-1 do begin
                    m = !values.f_nan
                    idz = where(z ge zbins[ib] and z lt zbins[ib+1], ngal)

                    if ngal ne 0 then begin
                        if p eq 0 then begin
                            m = mean((zphot.z_spec[idz] - z[idz])/(1.0 + zphot.z_spec[idz]))
                        endif else if p eq 1 then begin
                            m = 1.48*mad((zphot.z_spec[idz] - z[idz])/(1.0 + zphot.z_spec[idz]))
                        endif else if p eq 2 then begin
                            m = fraction_of((zphot.z_spec[idz] - z[idz])/(1.0 + zphot.z_spec[idz]) gt 0.15)*100
                        endif else if p eq 3 or p eq 4 then begin
                            if p eq 3 then dzth = 0.05
                            if p eq 4 then dzth = 0.15

                            ; Do the p(z) stacking in C++ because IDL is too slow
                            tmp_in = '/tmp/stack_list.fits'
                            tmp_out = out_stack_dir+'stack_'+zshow[iz]+'_b'+strn(ib)+'.fits'

                            if ~file_test(tmp_out) then begin
                                mwrfits, /create, {pz_file:result_dir+'pz/eazy_'+strna(idz+1)+'.pz', z_spec:zphot.z_spec[idz]}, tmp_in
                                spawn, './stack_pz '+tmp_in+' '+tmp_out
                            endif

                            stack = mrdfits(tmp_out, 1, /silent)

                            m = total(stack.pstack[where(abs(stack.dzgrid) lt dzth)])*100
                        endif
                    endif

                    my[2*ib+0] = m
                    my[2*ib+1] = m
                endfor

                idg = where(finite(my), cnt)
                if cnt ne 0 then begin
                    if p eq 0 then yr = [-max(abs(my[idg])) > (-0.5), +max(abs(my[idg])) < 0.5]
                    if p eq 1 then yr = [0, +max(abs(my[idg])) > 0.1 < 0.5]
                    if p eq 2 then yr = [0.0, max(my[idg]) > 20.0]
                    if p eq 3 then yr = [0.0, 100.0]
                    if p eq 4 then yr = [0.0, 100.0]
                endif else begin
                    if p eq 0 then yr = [-0.5, 0.5]
                    if p eq 1 then yr = [0, 0.5]
                    if p eq 2 then yr = [0.0, 20.0]
                    if p eq 3 then yr = [0.0, 100.0]
                    if p eq 4 then yr = [0.0, 100.0]
                endelse

                plot, /nodata, zx, my, xtit=xtit, xtickformat=xtickformat, ytit=ytit, $
                    pos=mplotpos(pos=p), yr=yr, xticklen=0.03, xr=[0,max(zbins)]

                if p eq 0 then oplotline, 0, 0, thick=3, line=0
                if p eq 0 then oplotline, 0, -0.002, thick=3, line=1
                if p eq 0 then oplotline, 0, +0.002, thick=3, line=1
                if p eq 1 then oplotline, 0, 0.05, thick=3, line=1
                if p eq 2 then oplotline, 0, 10, thick=3, line=1
                if p eq 3 then oplotline, 0, 68, thick=3, line=1
                if p eq 4 then oplotline, 0, 90, thick=3, line=1

                oplot, zx, my, thick=3, col='ff'x
            endfor
        endwhile
    endfor
endfor

end
