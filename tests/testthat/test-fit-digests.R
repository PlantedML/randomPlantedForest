# Small fits of every split structure, regression and classification, every classification loss, and purified fits
# (helper-fit-digests.R) must give these digests bitwise. They were computed with master at 6db527c and re-pinned
# where a fix moved the fits: the regression structures (`cur_trees_1`, and `cur_trees_2` / `res_trees` with
# `deterministic = TRUE`; #78), the purified fits (`purify(mode = 2)` on a leaf's lower bound, #76; `predict()`
# after a capped purification, #77), every classification fit (#84), the purified components of `purify(mode = 2)`
# as weighted rectangles, which add the same terms in another order (rounding level), and the three-class
# exponential fits (`loss3_exponential`, `pur_cls3_*`: SAMME, the multiclass exponential loss), the three-class
# L2 fit (`loss3_L2`: L2 classification scores a split with regression's arithmetic, which breaks a tie otherwise), and
# every logit and exponential fit (`cls_*`, `loss*_logit`, `loss*_exponential`, `pur_cls*` but their names: the
# losses' penalised step without clamps), and again every logit and exponential fit (the classification level,
# test-classif-level.R; with rpf(level = FALSE) they are master's). A change meant to leave
# every fit as it was keeps them all, and a change that moves fits re-pins exactly those.
#
# The digests depend on the platform's floating point and random integers (libm's exp and log draw the candidates),
# so they are checked only where they were computed: macOS on arm64. They do not depend on the BLAS / LAPACK that R
# uses (these fits call none; the level is computed without them, classif_level()), so R's reference BLAS (CRAN's
# macOS R), Accelerate and OpenBLAS give the same.

pinned_digests <- c(
  reg_leaves_det0_nt1 = "4d69f27123f4bf5b78bcb9228a33ae42",
  reg_leaves_det0_nt2 = "4d69f27123f4bf5b78bcb9228a33ae42",
  reg_leaves_det1_nt1 = "eb4a7b16933b5f1c88b2a6df1a8d5f63",
  cls_leaves = "bc74cd1927cf1a839d4ce242455bb3a9",
  reg_cur_trees_1_det0_nt1 = "f6bb739f4bc8010a420cda3570ca266d",
  reg_cur_trees_1_det0_nt2 = "f6bb739f4bc8010a420cda3570ca266d",
  reg_cur_trees_1_det1_nt1 = "f61d46f5d76d988905d0992f11b0be5c",
  cls_cur_trees_1 = "11719e5c548a12f4d1631048bac4a583",
  reg_cur_trees_2_det0_nt1 = "40d7845817eec7abf4547b51c2b25af7",
  reg_cur_trees_2_det0_nt2 = "40d7845817eec7abf4547b51c2b25af7",
  reg_cur_trees_2_det1_nt1 = "153661c73977f6b9a979e066f57fa3b0",
  cls_cur_trees_2 = "67c9e6f8d931b4da212e205f55272d35",
  reg_res_trees_det0_nt1 = "df541a66a06983f90c60bc21fe4cc317",
  reg_res_trees_det0_nt2 = "df541a66a06983f90c60bc21fe4cc317",
  reg_res_trees_det1_nt1 = "f3c531da7d664a9a8bff4ce0de0e1e9a",
  cls_res_trees = "0236fb00aedd61d08c4b41b74cd253a7",
  reg_hist_det0_nt1 = "bda13fd7237b1bf2bcb8358e75ad659e",
  reg_hist_det0_nt2 = "bda13fd7237b1bf2bcb8358e75ad659e",
  reg_hist_det1_nt1 = "42562c42ec0b77ec51ab7435eb1aefc8",
  cls_hist = "719991df65460831fa19ff46ac215275",
  loss2_L1 = "dba097acec118e611d815c3a827a9584",
  loss2_L2 = "85225648f3adc80ba30457e17b58bdcb",
  loss2_logit = "c56a5e11d60f2288d70c6586237aa3b0",
  loss2_exponential = "ebaa3818d1b9abba9b33c91b39faa9ce",
  loss3_L1 = "843dbb7d5d6bafd8eeea3bf3c67487d2",
  loss3_L2 = "c536e887470b870fe87d9304f9f55493",
  loss3_logit = "f6da6785a4c565d8371f8bfa9398fc3d",
  loss3_exponential = "12266abb766d52f722abbbc3bbeb7e59",
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
  pur_cls_components = "2dcb042bbbd344fabc8473fd487bbdf5",
  pur_cls_names = "794adb1c06fff2ed25b9afca6d9c654c",
  pur_cls_predict = "c8429429261b27535c675d9a4eb377ce",
  pur_cls3_components = "1a1bfe68d5da0fdac396ce5e54c3cfa3",
  pur_cls3_names = "28a2f21ce2d5f51e5a2bb15a2eee3069",
  pur_cls3_predict = "64d1b3fb84cfffba12f77636c788645c"
)

# The fits with rpf(level = FALSE): those of the version before the level, whose pins differed in the logit and
# exponential fits only (the rest, regression and L1 / L2, have no level)
pinned_digests_level_off <- replace(
  pinned_digests,
  c(
    "cls_leaves",
    "cls_cur_trees_1",
    "cls_cur_trees_2",
    "cls_res_trees",
    "cls_hist",
    "loss2_logit",
    "loss2_exponential",
    "loss3_logit",
    "loss3_exponential",
    "pur_cls_components",
    "pur_cls_predict",
    "pur_cls3_components",
    "pur_cls3_predict"
  ),
  c(
    "ea00222afc71241ba4cd8d68b1386b4e",
    "8342865d55b2727867ed04041ae31cc4",
    "41dd6f683e1cc359fb7120de1bb604f6",
    "cd41fe2592fa698efb964c42e0585ccd",
    "340c32f55abb8a82b31e86f60ce16b6b",
    "8c8fbbb8303640dc4c427d88cd5c6485",
    "8263e3e7d9a5db55c64b9b3c4e5f7f92",
    "83019105a5aa663469cebfba59029d89",
    "b1b1c73482efca7888f794065fee677b",
    "0cf62d6b6e8cadfae3d42cba3e5fb765",
    "14baa1eda3d55e5a1563d44b1537c377",
    "633193af37d38ee9f17b0d5df2cb3c23",
    "fef1a31c437d2efebbb14fdf3afc5866"
  )
)

expect_digests <- function(got, expected) {
  expect_identical(names(got), names(expected))
  for (nm in names(expected)) {
    expect_identical(got[[nm]], expected[[nm]], label = nm)
  }
}

skip_unless_pinned_platform <- function() {
  skip_on_cran()
  skip_if_not(
    Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64",
    "the pinned digests were computed on macOS arm64"
  )
}

test_that("small fits are bitwise as pinned", {
  skip_unless_pinned_platform()
  expect_digests(port_fit_digests(), pinned_digests)
})

test_that("with level = FALSE the small fits are bitwise those of the version before the level", {
  skip_unless_pinned_platform()
  expect_digests(port_fit_digests(port_fit_numbers_all(level = FALSE)), pinned_digests_level_off)
})
