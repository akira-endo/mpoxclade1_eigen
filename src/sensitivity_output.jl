if !(@isdefined bayesian) bayesian=true end
if !(@isdefined bcond) bcond = [0.1,0.01] end
if !(@isdefined estkeys) estkeys = [Symbol.(["1_" "2_"],"addmat",4:7)|>vec;:addmat_v;:addmat_w] end
# fit
kamituga2024_fit = output_sexual_fit(
    drc_kamituga,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_kamituga),
    estkeys = estkeys,
    bayesian=bayesian,
    bcond=bcond
  );

kivu2024_fit = output_sexual_fit(
    drc_kivu,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_kivu),
    estkeys = estkeys,
    bayesian=bayesian,
    bcond=bcond
    );

otherhz2024_fit = output_sexual_fit(
    drc_otherhz,
    zmb_sexual_skeleton = zmb2024_sexual,
    zmb_ref=last.(zmb2015_24_fit|>collect),
    drc_sexual_skeleton = drc2024_sexual,
    drc_ref=last.(drc2015_24_fit|>collect),
    dataplots = collapseplot(drc_otherhz),
    estkeys = estkeys,
    bayesian=bayesian,
    bcond=bcond
    );

burundi2024_fit = output_sexual_fit(
    burundi,
    zmb_sexual_skeleton = zmb2024_sexual_b,
    zmb_ref=zmb_ref,
    drc_sexual_skeleton = drc2024_sexual_b,
    drc_ref=drc_ref,
    dataplots = collapseplot(burundi),
    estkeys = estkeys,
    bayesian=bayesian,
    bcond=bcond
    );

if bayesian
    eig_kamituga=b_eigenanalysis(kamituga2024_fit.zmb_fit)
    eig_kivu=b_eigenanalysis(kivu2024_fit.zmb_fit)
    eig_otherhz=b_eigenanalysis(otherhz2024_fit.zmb_fit)
    eig_burundi=b_eigenanalysis(burundi2024_fit.zmb_fit)

    # fraction of R0 attributable to sexual transmisssion
    frac_R0_samples =[ (broadcast.(-,1, broadcast.(/,eig.eigval0, eig.eigval))) for eig in [eig_kivu, eig_kamituga, eig_otherhz, eig_burundi]]
    @show fracR0out = [median.(loc) for loc in frac_R0_samples]
else
    eig_kamituga=eigenanalysis(kamituga2024_fit.zmb_fit)
    eig_kivu=eigenanalysis(kivu2024_fit.zmb_fit)
    eig_otherhz=eigenanalysis(otherhz2024_fit.zmb_fit)
    eig_burundi=eigenanalysis(burundi2024_fit.zmb_fit);

    @show frac_R0out =[ (1 .-eig.eigval0./eig.eigval) for eig in [eig_kivu, eig_kamituga, eig_otherhz, eig_burundi]]
end
    