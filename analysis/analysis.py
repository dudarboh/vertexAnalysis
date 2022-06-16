import ROOT
import numpy as np

ROOT.gStyle.SetNdivisions(512)
ROOT.gStyle.SetNdivisions(512, "Y")
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)
ROOT.gStyle.SetOptStat("neMRou")
ROOT.gStyle.SetPaintTextFormat(".1f")
ROOT.gStyle.SetOptTitle(1)

ROOT.EnableImplicitMT(8)

#
# gInterpreter.Declare('''
#
#
#
#
# ''')

def plot_matrix():
    canvas = ROOT.TCanvas()
    # filename = "before_refit.root"
    filename = "after_refit.root"
    df = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/"+filename).Filter("has_kaon || has_proton")\
             .Define("n_true_tracks", "n_tracks - n_confused_tracks + n_missed_tracks")


    matrix = np.zeros((10, 10))

    data = df.AsNumpy(["n_tracks", "n_true_tracks"])

    for i, (n_tracks, n_true_tracks) in enumerate( zip(data["n_tracks"], data["n_true_tracks"]) ):
        if (0 < n_tracks < 10) and (0 < n_true_tracks < 10):
            matrix[n_tracks, n_true_tracks] += 1

    # Get Full matrix into histo for the cross check
    h = ROOT.TH2F("h", ";N_{reco} tracks;N_{true} tracks", 8, 2, 10, 8, 2, 10)
    h.SetStats(0)
    ROOT.gStyle.SetPaintTextFormat(".2f")

    for i in range(8):
        for j in range(8):
            h.SetBinContent(i+1, j+1, matrix[i+2, j+2])

    for r in range(1, 9):
        norm = sum( [h.GetBinContent(i, r) for i in range(1, 11)] )
        if norm == 0:
            continue
        for c in range(1, 9):
            h.SetBinContent(c, r, h.GetBinContent(c, r) / norm)
    h.Draw("colz text")
    canvas.Update()
    input("wait")






## BEFORE REFIT

canvas = ROOT.TCanvas("dtoip_kaons_3gev_norefit_cut")
canvas.SetGridx()
canvas.SetGridy()
canvas.Divide(2, 2)

canvas.cd(1)
df1 = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/before_refit_with_cut.root")
# df1 = df1.Filter("n_true_kaons == 1 && p_true_kaon.r() < 3.").Define("dtoip", "(pos_true - ip_pos).r()")
df1 = df1.Define("dtoip", "(pos_true - ip_pos).r()")
h_missed1 = df1.Filter("is_missed == 1 ").Histo1D(("h_missed1", "Missed; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_missed1.SetFillColor(ROOT.kRed+1)
h_confused1 = df1.Filter("is_confused == 1  && n_confused_tracks != 0 ").Histo1D(("h_confused1", "Confused; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_confused1.SetFillColor(ROOT.kGreen+2)
h_matched1 = df1.Filter("is_matched == 1  || (is_confused == 1 && n_confused_tracks == 0)").Histo1D(("h_matched1", "Matched; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_matched1.SetFillColor(ROOT.kBlue+1)
hs1 = ROOT.THStack()
hs1.SetTitle("kaons 3 GeV before_refit_with_cut;d to IP (mm); N vertices")
hs1.Add(h_matched1.GetPtr())
hs1.Add(h_confused1.GetPtr())
hs1.Add(h_missed1.GetPtr())
hs1.Draw()
#hs1.SetMaximum(6000)
ROOT.gPad.BuildLegend()
# input("wait")

canvas.cd(2)
df2 = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/after_refit_with_cut.root")
# df2 = df2.Filter("n_true_kaons == 1 && p_true_kaon.r() < 3.").Define("dtoip", "(pos_true - ip_pos).r()")
df2 = df2.Define("dtoip", "(pos_true - ip_pos).r()")
h_missed2 = df2.Filter("is_missed == 1 ").Histo1D(("h_missed2", "Missed; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_missed2.SetFillColor(ROOT.kRed+1)
h_confused2 = df2.Filter("is_confused == 1  && n_confused_tracks != 0 ").Histo1D(("h_confused2", "Confused; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_confused2.SetFillColor(ROOT.kGreen+2)
h_matched2 = df2.Filter("is_matched == 1  || (is_confused == 1 && n_confused_tracks == 0)").Histo1D(("h_matched2", "Matched; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_matched2.SetFillColor(ROOT.kBlue+1)
hs2 = ROOT.THStack()
hs2.SetTitle("kaons 3 GeV after_refit_with_cut;d to IP (mm); N vertices")
hs2.Add(h_matched2.GetPtr())
hs2.Add(h_confused2.GetPtr())
hs2.Add(h_missed2.GetPtr())
hs2.Draw()
#hs2.SetMaximum(6000)
# input("wait")

canvas.cd(3)
df3 = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/before_refit_without_cut.root")
# df3 = df3.Filter("n_true_kaons == 1 && p_true_kaon.r() < 3.").Define("dtoip", "(pos_true - ip_pos).r()")
df3 = df3.Define("dtoip", "(pos_true - ip_pos).r()")
h_missed3 = df3.Filter("is_missed == 1 ").Histo1D(("h_missed3", "Missed; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_missed3.SetFillColor(ROOT.kRed+1)
h_confused3 = df3.Filter("is_confused == 1  && n_confused_tracks != 0 ").Histo1D(("h_confused3", "Confused; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_confused3.SetFillColor(ROOT.kGreen+2)
h_matched3 = df3.Filter("is_matched == 1  || (is_confused == 1 && n_confused_tracks == 0)").Histo1D(("h_matched3", "Matched; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_matched3.SetFillColor(ROOT.kBlue+1)
hs3 = ROOT.THStack()
hs3.SetTitle("kaons 3 GeV before_refit_without_cut;d to IP (mm); N vertices")
hs3.Add(h_matched3.GetPtr())
hs3.Add(h_confused3.GetPtr())
hs3.Add(h_missed3.GetPtr())
hs3.Draw()
#hs3.SetMaximum(6000)
# input("wait")

canvas.cd(4)
df4 = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/after_refit_without_cut.root")
# df4 = df4.Filter("n_true_kaons == 1 && p_true_kaon.r() < 3.").Define("dtoip", "(pos_true - ip_pos).r()")
df4 = df4.Define("dtoip", "(pos_true - ip_pos).r()")
h_missed4 = df4.Filter("is_missed == 1 ").Histo1D(("h_missed4", "Missed; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_missed4.SetFillColor(ROOT.kRed+1)
h_confused4 = df4.Filter("is_confused == 1  && n_confused_tracks != 0 ").Histo1D(("h_confused4", "Confused; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_confused4.SetFillColor(ROOT.kGreen+2)
h_matched4 = df4.Filter("is_matched == 1  || (is_confused == 1 && n_confused_tracks == 0)").Histo1D(("h_matched4", "Matched; |#vec{r}_{true} - #vec{r}_{IP}| (mm); N vertices", 1000, -1., 10.), "dtoip")
h_matched4.SetFillColor(ROOT.kBlue+1)
hs4 = ROOT.THStack()
hs4.SetTitle("kaons 3 GeV after_refit_without_cut;d to IP (mm); N vertices")
hs4.Add(h_matched4.GetPtr())
hs4.Add(h_confused4.GetPtr())
hs4.Add(h_missed4.GetPtr())
hs4.Draw()
#hs4.SetMaximum(6000)
canvas.Update()
input("wait")

# print("BEFORE REFIT NO CUT: N matched ", n_matched_old.GetValue(), "(", 100.*n_matched_old.GetValue() / n_total_old.GetValue() , ")")
# print("BEFORE REFIT NO CUT: N missed ", n_missed_old.GetValue(), "(", 100.*n_missed_old.GetValue() / n_total_old.GetValue() , ")")
# print("BEFORE REFIT NO CUT: N confused ", n_confused_old.GetValue(), "(", 100.*n_confused_old.GetValue() / n_total_old.GetValue() , ")")
# 
# print("AFTER REFIT: N matched ", n_matched_new.GetValue(), "(", 100.*n_matched_new.GetValue() / n_total_new.GetValue() , ")")
# print("AFTER REFIT: N missed ", n_missed_new.GetValue(), "(", 100.*n_missed_new.GetValue() / n_total_new.GetValue() , ")")
# print("AFTER REFIT: N confused ", n_confused_new.GetValue(), "(", 100.*n_confused_new.GetValue() / n_total_new.GetValue() , ")")




# h_ntracks = df_old.Histo1D(("h_ntracks", "title;n tracks;y", 10, 0, 10), "n_true_tracks")
# h_ntracks.Draw()
# input("wait")
# 

n1_total = df1.Count()
n1_matched = df1.Filter("is_matched == 1 || (is_confused == 1 && n_confused_tracks == 0)").Count()
n1_missed = df1.Filter("is_missed == 1").Count()
n1_confused = df1.Filter("is_confused == 1 && n_confused_tracks != 0").Count()

# 
# 
# 
# 
# 
# df_new = ROOT.RDataFrame("VertexAnalysis", "/nfs/dust/ilc/user/dudarboh/tof/refit/after_refit_without_cut.root")
# df_new = df_new.Define("mom", "p_true_proton.r()").Filter("n_true_kaons == 1 && mom < 3.")
# n_total_new = df_new.Count()
# n_matched_new = df_new.Filter("is_matched == 1").Count()
# n_missed_new = df_new.Filter("is_missed == 1").Count()
# n_confused_new = df_new.Filter("is_confused == 1").Count()




# df_new = df_new.Define("mom", "p_reco_proton.r()").Filter(" n_reco_protons > 0 && mom < 1.5")
# df_new = df_new.Define("mom", "p_reco_proton.r()").Filter(" n_reco_protons > 0 && mom < 1.5")
# df_new = df_new.Define("mom", "p_reco_proton.r()").Filter("n_reco_protons > 0")

################################### Position diff
# h_old_dtoip = df_old.Define("diff_r", "(pos_true - ip_pos).r()*1000").Histo1D(("h_old_r", "w/o refit; |#vec{r}_{reco} - #vec{r}_{MC}| (#mum); N vertices", 100, -1., 200.), "diff_r")
# h_new_dtoip = df_new.Define("diff_r", "(pos_true - ip_pos).r()*1000").Histo1D(("h_new_r", "w/ refit; |#vec{r}_{reco} - #vec{r}_{MC}| (#mum); N vertices", 100, -1., 200.), "diff_r")
# 
# 
# df_old = df_old.Filter("is_matched == 1")
# df_new = df_new.Filter("is_matched == 1")
# h_old_x = df_old.Define("diff_x", "(pos_reco - pos_true).x()*1000").Histo1D(("h_old_x", "w/o refit; x_{reco} - x_{MC} (#mum); N vertices", 100, -200., 200.), "diff_x")
# h_new_x = df_new.Define("diff_x", "(pos_reco - pos_true).x()*1000").Histo1D(("h_new_x", "w/ refit; x_{reco} - x_{MC} (#mum); N vertices", 100, -200., 200.), "diff_x")
# h_old_y = df_old.Define("diff_y", "(pos_reco - pos_true).y()*1000").Histo1D(("h_old_y", "w/o refit; y_{reco} - y_{MC} (#mum); N vertices", 100, -200., 200.), "diff_y")
# h_new_y = df_new.Define("diff_y", "(pos_reco - pos_true).y()*1000").Histo1D(("h_new_y", "w/ refit; y_{reco} - y_{MC} (#mum); N vertices", 100, -200., 200.), "diff_y")
# h_old_z = df_old.Define("diff_z", "(pos_reco - pos_true).z()*1000").Histo1D(("h_old_z", "w/o refit; z_{reco} - z_{MC} (#mum); N vertices", 100, -200., 200.), "diff_z")
# h_new_z = df_new.Define("diff_z", "(pos_reco - pos_true).z()*1000").Histo1D(("h_new_z", "w/ refit; z_{reco} - z_{MC} (#mum); N vertices", 100, -200., 200.), "diff_z")
# h_old_r = df_old.Define("diff_r", "(pos_reco - pos_true).r()*1000").Histo1D(("h_old_r", "w/o refit; |#vec{r}_{reco} - #vec{r}_{MC}| (#mum); N vertices", 100, -1., 200.), "diff_r")
# h_new_r = df_new.Define("diff_r", "(pos_reco - pos_true).r()*1000").Histo1D(("h_new_r", "w/ refit; |#vec{r}_{reco} - #vec{r}_{MC}| (#mum); N vertices", 100, -1., 200.), "diff_r")



# h_old_x = df_old.Define("diff_x", "(pos_reco - pos_true).x()/sigma_x").Histo1D(("h_old_x", "w/o refit; (x_{reco} - x_{MC}) / #sigma_{x} ; N vertices", 100, -10, 10), "diff_x")
# h_new_x = df_new.Define("diff_x", "(pos_reco - pos_true).x()/sigma_x").Histo1D(("h_new_x", "w/ refit; (x_{reco} - x_{MC}) /#sigma_{x} ; N vertices", 100, -10, 10), "diff_x")
# h_old_y = df_old.Define("diff_y", "(pos_reco - pos_true).y()/sigma_y").Histo1D(("h_old_y", "w/o refit; (y_{reco} - y_{MC}) /#sigma_{y} ; N vertices", 100, -10, 10), "diff_y")
# h_new_y = df_new.Define("diff_y", "(pos_reco - pos_true).y()/sigma_y").Histo1D(("h_new_y", "w/ refit; (y_{reco} - y_{MC}) /#sigma_{y} ; N vertices", 100, -10, 10), "diff_y")
# h_old_z = df_old.Define("diff_z", "(pos_reco - pos_true).z()/sigma_z").Histo1D(("h_old_z", "w/o refit; (z_{reco} - z_{MC}) /#sigma_{z} ; N vertices", 100, -10, 10), "diff_z")
# h_new_z = df_new.Define("diff_z", "(pos_reco - pos_true).z()/sigma_z").Histo1D(("h_new_z", "w/ refit; (z_{reco} - z_{MC}) /#sigma_{z} ; N vertices", 100, -10, 10), "diff_z")
# h_old_r = df_old.Define("diff_r", "(pos_reco - pos_true).r()/sqrt(sigma_x*sigma_x + sigma_y*sigma_y + sigma_z*sigma_z)").Histo1D(("h_old_r", "w/o refit; |#vec{r}_{reco} - #vec{r}_{MC}| / #sigma_{r} ; N vertices", 100, 0, 10), "diff_r")
# h_new_r = df_new.Define("diff_r", "(pos_reco - pos_true).r()/sqrt(sigma_x*sigma_x + sigma_y*sigma_y + sigma_z*sigma_z)").Histo1D(("h_new_r", "w/ refit; |#vec{r}_{reco} - #vec{r}_{MC}| / #sigma_{r} ; N vertices", 100, 0, 10), "diff_r")

# print("BEFORE REFIT: N matched ", n_matched_old.GetValue(), "(", 100.*n_matched_old.GetValue() / n_total_old.GetValue() , ")")
# print("BEFORE REFIT: N missed ", n_missed_old.GetValue(), "(", 100.*n_missed_old.GetValue() / n_total_old.GetValue() , ")")
# print("BEFORE REFIT: N confused ", n_confused_old.GetValue(), "(", 100.*n_confused_old.GetValue() / n_total_old.GetValue() , ")")
# 
# print("AFTER REFIT: N matched ", n_matched_new.GetValue(), "(", 100.*n_matched_new.GetValue() / n_total_new.GetValue() , ")")
# print("AFTER REFIT: N missed ", n_missed_new.GetValue(), "(", 100.*n_missed_new.GetValue() / n_total_new.GetValue() , ")")
# print("AFTER REFIT: N confused ", n_confused_new.GetValue(), "(", 100.*n_confused_new.GetValue() / n_total_new.GetValue() , ")")


# # h_old_x.Scale(1. / n_old_vertices.GetValue())
# # h_new_x.Scale(1. / n_new_vertices.GetValue())
# h_old_x.Draw("histo")
# h_new_x.Draw("histo sames")
# h_old_x.SetMaximum(1.05*max(h_old_x.GetMaximum(), h_new_x.GetMaximum()))
# # print(h_old_x.FindObject("stats"))
# # stats_old_x = h_old_x.FindObject("stats")
# # stats_old_x.SetX1NDC(0.13)
# # stats_old_x.SetX2NDC(0.28)
# h_new_x.SetLineColor(2)
# # stats_new_x = h_new_x.FindObject("stats")
# # stats_new_x.SetLineColor(2)
# # stats_new_x.SetTextColor(2)
# 
# 
# canvas.cd(2)
# # h_old_y.Scale(1. / n_old_vertices.GetValue())
# # h_new_y.Scale(1. / n_new_vertices.GetValue())
# h_old_y.Draw("histo")
# h_new_y.Draw("histo sames")
# h_old_y.SetMaximum(1.05*max(h_old_y.GetMaximum(), h_new_y.GetMaximum()))
# h_new_y.SetLineColor(2)
# 
# canvas.cd(3)
# # h_old_z.Scale(1. / n_old_vertices.GetValue())
# # h_new_z.Scale(1. / n_new_vertices.GetValue())
# h_old_z.Draw("histo")
# h_new_z.Draw("histo sames")
# h_old_z.SetMaximum(1.05*max(h_old_z.GetMaximum(), h_new_z.GetMaximum()))
# h_new_z.SetLineColor(2)
# 
# canvas.cd(4)
# # h_old_r.Scale(1. / n_old_vertices.GetValue())
# # h_new_r.Scale(1. / n_new_vertices.GetValue())
# h_old_r.Draw("histo")
# h_new_r.Draw("histo sames")
# h_old_r.SetMaximum(1.05*max(h_old_r.GetMaximum(), h_new_r.GetMaximum()))
# h_new_r.SetLineColor(2)
# 
# 
# h_old_x.SetTitle("X residuals")
# h_old_y.SetTitle("Y residuals")
# h_old_z.SetTitle("Z residuals")
# h_old_r.SetTitle("R residuals")


######################################################## PUll vs momentum

# h_old_x = df_old.Define("diff_x", "(pos_reco - pos_true).x()/sigma_x").Histo2D(("h_old_x", "w/o refit;  p (GeV);(x_{reco} - x_{MC}) / #sigma_{x}", 20, 0, 5, 100, -10, 10), "mom", "diff_x")
# h_new_x = df_new.Define("diff_x", "(pos_reco - pos_true).x()/sigma_x").Histo2D(("h_new_x", "w/ refit; p (GeV); (x_{reco} - x_{MC}) /#sigma_{x}", 20, 0, 5, 100, -10, 10), "mom", "diff_x")
# h_old_y = df_old.Define("diff_y", "(pos_reco - pos_true).y()/sigma_y").Histo2D(("h_old_y", "w/o refit; p (GeV); (y_{reco} - y_{MC}) /#sigma_{y}", 20, 0, 5, 100, -10, 10), "mom", "diff_y")
# h_new_y = df_new.Define("diff_y", "(pos_reco - pos_true).y()/sigma_y").Histo2D(("h_new_y", "w/ refit; p (GeV); (y_{reco} - y_{MC}) /#sigma_{y}", 20, 0, 5, 100, -10, 10), "mom", "diff_y")
# h_old_z = df_old.Define("diff_z", "(pos_reco - pos_true).z()/sigma_z").Histo2D(("h_old_z", "w/o refit; p (GeV); (z_{reco} - z_{MC}) /#sigma_{z}", 20, 0, 5, 100, -10, 10), "mom", "diff_z")
# h_new_z = df_new.Define("diff_z", "(pos_reco - pos_true).z()/sigma_z").Histo2D(("h_new_z", "w/ refit; p (GeV); (z_{reco} - z_{MC}) /#sigma_{z}", 20, 0, 5, 100, -10, 10), "mom", "diff_z")

# h_old_x.Draw("colz")
# canvas2 = ROOT.TCanvas()
# canva2.cd()
# prof = h_old_z.ProfileX()
# prof.BuildOptions(-10, 10, "s")
# prof.Draw()
# prof_new = h_new_z.ProfileX()
# prof_new.BuildOptions(-10, 10, "s")
# prof_new.Draw("same")
# prof_new.SetLineColor(2)
################################### N confused tracks
# h_old = df_old.Histo1D(("h_old", "w/o refit; N confused tracks; N build up vertices (%)", 10, 0, 10), "n_confused_tracks")
# h_new = df_new.Histo1D(("h_new", "w/ refit; N confused tracks; N build up vertices (%)", 10, 0, 10), "n_confused_tracks")
# h_old.Scale(1./h_old.GetEntries())
# h_new.Scale(1./h_new.GetEntries())
# h_old.Draw("histo")
# h_new.Draw("histo sames")
# h_new.SetLineColor(2)

################################### 2D confusion plot
# plot_matrix()

################################### parent hadrons
# h_old = df_old.Histo1D(("h_old", "w/o refit; p (GeV); N vertices", 40, 0, 20), "mom")
# h_new = df_new.Histo1D(("h_new", "w/ refit; p (GeV); N vertices", 40, 0, 20), "mom")
# h_old.Scale(1./n_old_vertices.GetValue())
# h_new.Scale(1./n_new_vertices.GetValue())
# print("N vertices OLD: ", n_old_vertices.GetValue())
# print("N vertices NEW: ", n_new_vertices.GetValue())
# h_old.Draw("histo")
# h_new.Draw("histo sames")
# h_new.SetLineColor(2)




################################### Position diff sigma = pulls
# h_old_x = df_old.Define("diff_x", "(reco_vtx_pos - mc_vtx_pos).x()/std::sqrt(cov_xx)").Histo1D(("h_old_x", "w/o refit; x_{reco} - x_{MC} / #sigma_{xx}; N build up vertices", 400, -10., 10.), "diff_x")
# h_new_x = df_new.Define("diff_x", "(reco_vtx_pos - mc_vtx_pos).x()/std::sqrt(cov_xx)").Histo1D(("h_new_x", "w/ refit; x_{reco} - x_{MC} / #sigma_{xx}; N build up vertices", 400, -10., 10.), "diff_x")
# h_old_y = df_old.Define("diff_y", "(reco_vtx_pos - mc_vtx_pos).y()/std::sqrt(cov_yy)").Histo1D(("h_old_y", "w/o refit; y_{reco} - y_{MC} / #sigma_{yy}; N build up vertices", 400, -10., 10.), "diff_y")
# h_new_y = df_new.Define("diff_y", "(reco_vtx_pos - mc_vtx_pos).y()/std::sqrt(cov_yy)").Histo1D(("h_new_y", "w/ refit; y_{reco} - y_{MC} / #sigma_{yy}; N build up vertices", 400, -10., 10.), "diff_y")
# h_old_z = df_old.Define("diff_z", "(reco_vtx_pos - mc_vtx_pos).z()/std::sqrt(cov_zz)").Histo1D(("h_old_z", "w/o refit; z_{reco} - z_{MC} / #sigma_{zz}; N build up vertices", 400, -10., 10.), "diff_z")
# h_new_z = df_new.Define("diff_z", "(reco_vtx_pos - mc_vtx_pos).z()/std::sqrt(cov_zz)").Histo1D(("h_new_z", "w/ refit; z_{reco} - z_{MC} / #sigma_{zz}; N build up vertices", 400, -10., 10.), "diff_z")
# h_old_r = df_old.Define("diff_r", "(reco_vtx_pos - mc_vtx_pos).r()/std::sqrt(cov_xx+cov_yy+cov_zz)").Histo1D(("h_old_r", "w/o refit; |#vec{r}_{reco} - #vec{r}_{MC}| /#sigma_{tot}; N build up vertices", 400, -10., 10.), "diff_r")
# h_new_r = df_new.Define("diff_r", "(reco_vtx_pos - mc_vtx_pos).r()/std::sqrt(cov_xx+cov_yy+cov_zz)").Histo1D(("h_new_r", "w/ refit; |#vec{r}_{reco} - #vec{r}_{MC}| /#sigma_{tot}; N build up vertices", 400, -10., 10.), "diff_r")
# canvas.Divide(2, 2)
# 
# canvas.cd(1)
# h_old_x.Draw()
# h_new_x.Draw("sames")
# h_new_x.SetLineColor(2)
# 
# canvas.cd(2)
# h_old_y.Draw()
# h_new_y.Draw("sames")
# h_new_y.SetLineColor(2)
# 
# canvas.cd(3)
# h_old_z.Draw()
# h_new_z.Draw("sames")
# h_new_z.SetLineColor(2)
# 
# canvas.cd(4)
# # h_old_r.Draw()
# # h_new_r.Draw("sames")
# # h_new_r.SetLineColor(2)
# 
# 
# h_old_x.SetTitle("X normalized residual")
# h_old_y.SetTitle("Y normalized residual")
# h_old_z.SetTitle("Z normalized residual")
# h_old_r.SetTitle("R normalized residual")





# canvas.BuildLegend()
canvas.Update()

input("wait")
