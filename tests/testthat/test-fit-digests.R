# Small fits of every split structure, regression and classification, every classification loss, and purified fits
# (helper-fit-digests.R) must give these digests bitwise. They were computed with master at 6db527c and re-pinned
# where a fix moved the fits: the regression structures (`cur_trees_1`, and `cur_trees_2` / `res_trees` with
# `deterministic = TRUE`; #78), the purified fits (`purify(mode = 2)` on a leaf's lower bound, #76; `predict()`
# after a capped purification, #77) and every classification fit (#84). A change meant to leave every fit as it was
# keeps them all, and a change that moves fits re-pins exactly those.
#
# The digests depend on the platform's floating point and random integers (libm's exp and log draw the candidates),
# so they are checked only where they were computed: macOS on arm64.

test_that("small fits are bitwise as pinned", {
  skip_on_cran()
  skip_if_not(
    Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64",
    "the pinned digests were computed on macOS arm64"
  )
  expected <- c(
    reg_leaves_det0_nt1 = "4d69f27123f4bf5b78bcb9228a33ae42",
    reg_leaves_det0_nt2 = "4d69f27123f4bf5b78bcb9228a33ae42",
    reg_leaves_det1_nt1 = "eb4a7b16933b5f1c88b2a6df1a8d5f63",
    cls_leaves = "4276783f4e5b1fdf4107d0dfcbe36e07",
    reg_cur_trees_1_det0_nt1 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det0_nt2 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det1_nt1 = "f61d46f5d76d988905d0992f11b0be5c",
    cls_cur_trees_1 = "93b6de0783486c8109dc3d1742a0db47",
    reg_cur_trees_2_det0_nt1 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det0_nt2 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det1_nt1 = "153661c73977f6b9a979e066f57fa3b0",
    cls_cur_trees_2 = "6aa6e96292eaee7c4c643f03bff629cc",
    reg_res_trees_det0_nt1 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det0_nt2 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det1_nt1 = "f3c531da7d664a9a8bff4ce0de0e1e9a",
    cls_res_trees = "8ed48ff55cf9f953fea3ee10c310151e",
    reg_hist_det0_nt1 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det0_nt2 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det1_nt1 = "42562c42ec0b77ec51ab7435eb1aefc8",
    cls_hist = "7bfec596388d019976ed704ddf08ffc2",
    loss2_L1 = "dba097acec118e611d815c3a827a9584",
    loss2_L2 = "85225648f3adc80ba30457e17b58bdcb",
    loss2_logit = "fff9199d1334d380b87491678104fa2a",
    loss2_exponential = "de10362d03c0d12ab8a5fd66e2d62d86",
    loss3_L1 = "843dbb7d5d6bafd8eeea3bf3c67487d2",
    loss3_L2 = "05d608d4ff10164f09ff56353c399223",
    loss3_logit = "305b7853bd96e566572be447023b34d1",
    loss3_exponential = "04b5ee5e76a1ae620126dd83c68a3fcc",
    pur_reg_components = "78dc819b51ff0c8dc8fe34fc481ef773",
    pur_reg_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_predict = "b5fbdc9b147dc349dc6126d46b9b0ddb",
    pur_reg_nt2_components = "78dc819b51ff0c8dc8fe34fc481ef773",
    pur_reg_nt2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_nt2_predict = "b5fbdc9b147dc349dc6126d46b9b0ddb",
    pur_reg_mode1_components = "30432c23c98b6558ef55e04879542d36",
    pur_reg_mode1_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_mode1_predict = "481c2302637dffaa23e872947fb6fedc",
    pur_reg_cap2_components = "597ebd88a997ceb66a92bf132deb0649",
    pur_reg_cap2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_cap2_predict = "0cc253088d2a69994515aeead292a80c",
    pur_cls_components = "6553d4cc52b83c34210ce37dbecf07e0",
    pur_cls_names = "794adb1c06fff2ed25b9afca6d9c654c",
    pur_cls_predict = "d78fc63f666c04dc4c01f46e5032a0a7",
    pur_cls3_components = "c384030fe2ad549ca38563cf20165714",
    pur_cls3_names = "28a2f21ce2d5f51e5a2bb15a2eee3069",
    pur_cls3_predict = "d846be30b5b61e1b24f21ffcb43486b7"
  )
  got <- port_fit_digests()
  expect_identical(names(got), names(expected))
  for (nm in names(expected)) {
    expect_identical(got[[nm]], expected[[nm]], label = nm)
  }
})
