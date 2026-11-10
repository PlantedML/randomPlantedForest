# Small fits of every split structure, regression and classification, every classification loss, and purified fits
# (helper-fit-digests.R) must give these digests bitwise. They were computed with master at 6db527c and re-pinned
# where a fix moved the fits: the regression structures (`cur_trees_1`, and `cur_trees_2` / `res_trees` with
# `deterministic = TRUE`; #78) and the purified fits (`purify(mode = 2)` on a leaf's lower bound, #76; `predict()`
# after a capped purification, #77). A change meant to leave every fit as it was keeps them all, and a change that
# moves fits re-pins exactly those.
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
    cls_leaves = "6bf49f16e4383bef971ea7ea2e09190f",
    reg_cur_trees_1_det0_nt1 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det0_nt2 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det1_nt1 = "f61d46f5d76d988905d0992f11b0be5c",
    cls_cur_trees_1 = "aeba1c4138772364f54675f736f6e03a",
    reg_cur_trees_2_det0_nt1 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det0_nt2 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det1_nt1 = "153661c73977f6b9a979e066f57fa3b0",
    cls_cur_trees_2 = "aeba1c4138772364f54675f736f6e03a",
    reg_res_trees_det0_nt1 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det0_nt2 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det1_nt1 = "f3c531da7d664a9a8bff4ce0de0e1e9a",
    cls_res_trees = "8196f89f59b1442cc5e04af2086a0c74",
    reg_hist_det0_nt1 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det0_nt2 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det1_nt1 = "42562c42ec0b77ec51ab7435eb1aefc8",
    cls_hist = "ef37a8e9056997842b82844e25c5c818",
    loss2_L1 = "ef90ecd989b6cc7f226a860a36f611a6",
    loss2_L2 = "d828ad1633549ada3e9088d3b127fa13",
    loss2_logit = "0f25e64f4978588ffdaf5b14baf68dc4",
    loss2_exponential = "a2936e221f72931b4912bce3aaf2cec1",
    loss3_L1 = "0c0fa24236f44db4ccad3cc4178549aa",
    loss3_L2 = "1252eba810b3193bc25ca20b4e163063",
    loss3_logit = "03bf37adc2da3412bbe3671bed4baade",
    loss3_exponential = "eb0b7f309df5f548a25f3903904b0e44",
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
    pur_cls_components = "b0ac80a5f8e08745874805da4ef6161d",
    pur_cls_names = "794adb1c06fff2ed25b9afca6d9c654c",
    pur_cls_predict = "31076f20aabc4afdf91325c800f83805",
    pur_cls3_components = "0581fdc8a24d4e0ed4caf18a13477cc6",
    pur_cls3_names = "28a2f21ce2d5f51e5a2bb15a2eee3069",
    pur_cls3_predict = "3879324b0e8e2fa0e15493fe2133ee73"
  )
  got <- port_fit_digests()
  expect_identical(names(got), names(expected))
  for (nm in names(expected)) {
    expect_identical(got[[nm]], expected[[nm]], label = nm)
  }
})
