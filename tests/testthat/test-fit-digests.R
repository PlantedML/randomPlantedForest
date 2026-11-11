# Small fits of every split structure, regression and classification, every classification loss, and purified fits
# (helper-fit-digests.R) must give these digests bitwise. They were computed with master at 6db527c and re-pinned
# where a fix moved the fits: the regression structures (`cur_trees_1`, and `cur_trees_2` / `res_trees` with
# `deterministic = TRUE`; #78), the purified fits (`purify(mode = 2)` on a leaf's lower bound, #76; `predict()`
# after a capped purification, #77), every classification fit (#84), the purified components of `purify(mode = 2)`
# as weighted rectangles, which add the same terms in another order (rounding level), and the three-class
# exponential fits (`loss3_exponential`, `pur_cls3_*`: SAMME, the multiclass exponential loss), the three-class
# L2 fit (`loss3_L2`: L2 classification scores a split with regression's arithmetic, which breaks a tie otherwise), and
# every logit and exponential fit (`cls_*`, `loss*_logit`, `loss*_exponential`, `pur_cls*` but their names: the
# losses' penalised step without clamps). A change meant to leave every fit as it was keeps them all, and a change
# that moves fits re-pins exactly those.
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
    cls_leaves = "ea00222afc71241ba4cd8d68b1386b4e",
    reg_cur_trees_1_det0_nt1 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det0_nt2 = "f6bb739f4bc8010a420cda3570ca266d",
    reg_cur_trees_1_det1_nt1 = "f61d46f5d76d988905d0992f11b0be5c",
    cls_cur_trees_1 = "8342865d55b2727867ed04041ae31cc4",
    reg_cur_trees_2_det0_nt1 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det0_nt2 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det1_nt1 = "153661c73977f6b9a979e066f57fa3b0",
    cls_cur_trees_2 = "41dd6f683e1cc359fb7120de1bb604f6",
    reg_res_trees_det0_nt1 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det0_nt2 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det1_nt1 = "f3c531da7d664a9a8bff4ce0de0e1e9a",
    cls_res_trees = "cd41fe2592fa698efb964c42e0585ccd",
    reg_hist_det0_nt1 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det0_nt2 = "bda13fd7237b1bf2bcb8358e75ad659e",
    reg_hist_det1_nt1 = "42562c42ec0b77ec51ab7435eb1aefc8",
    cls_hist = "340c32f55abb8a82b31e86f60ce16b6b",
    loss2_L1 = "dba097acec118e611d815c3a827a9584",
    loss2_L2 = "85225648f3adc80ba30457e17b58bdcb",
    loss2_logit = "8c8fbbb8303640dc4c427d88cd5c6485",
    loss2_exponential = "8263e3e7d9a5db55c64b9b3c4e5f7f92",
    loss3_L1 = "843dbb7d5d6bafd8eeea3bf3c67487d2",
    loss3_L2 = "c536e887470b870fe87d9304f9f55493",
    loss3_logit = "83019105a5aa663469cebfba59029d89",
    loss3_exponential = "b1b1c73482efca7888f794065fee677b",
    pur_reg_components = "8bf9507e8cb6ea716ec1315c58a5de81",
    pur_reg_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_predict = "e6f84aa30eedf710528f06e78974a9ea",
    pur_reg_nt2_components = "8bf9507e8cb6ea716ec1315c58a5de81",
    pur_reg_nt2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_nt2_predict = "e6f84aa30eedf710528f06e78974a9ea",
    pur_reg_mode1_components = "30432c23c98b6558ef55e04879542d36",
    pur_reg_mode1_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_mode1_predict = "481c2302637dffaa23e872947fb6fedc",
    pur_reg_cap2_components = "d45cc7f004f85b3d25ece7b214ee21d4",
    pur_reg_cap2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_cap2_predict = "0cc253088d2a69994515aeead292a80c",
    pur_cls_components = "0cf62d6b6e8cadfae3d42cba3e5fb765",
    pur_cls_names = "794adb1c06fff2ed25b9afca6d9c654c",
    pur_cls_predict = "14baa1eda3d55e5a1563d44b1537c377",
    pur_cls3_components = "633193af37d38ee9f17b0d5df2cb3c23",
    pur_cls3_names = "28a2f21ce2d5f51e5a2bb15a2eee3069",
    pur_cls3_predict = "fef1a31c437d2efebbb14fdf3afc5866"
  )
  got <- port_fit_digests()
  expect_identical(names(got), names(expected))
  for (nm in names(expected)) {
    expect_identical(got[[nm]], expected[[nm]], label = nm)
  }
})
