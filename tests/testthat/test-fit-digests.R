# Small fits of every split structure, regression and classification, every classification loss, and purified fits
# (helper-fit-digests.R) must give these digests bitwise. They were computed with master at 6db527c: a change meant
# to leave every fit as it was keeps them all, and a change that moves fits re-pins exactly those.
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
    reg_cur_trees_1_det0_nt1 = "d80f92b149f431d19a00445e489cf55a",
    reg_cur_trees_1_det0_nt2 = "d80f92b149f431d19a00445e489cf55a",
    reg_cur_trees_1_det1_nt1 = "6d67e66eeec3a043d4e9bd1a87d64ef0",
    cls_cur_trees_1 = "aeba1c4138772364f54675f736f6e03a",
    reg_cur_trees_2_det0_nt1 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det0_nt2 = "40d7845817eec7abf4547b51c2b25af7",
    reg_cur_trees_2_det1_nt1 = "bf8d498986839bbb56c18df11e8e6c14",
    cls_cur_trees_2 = "aeba1c4138772364f54675f736f6e03a",
    reg_res_trees_det0_nt1 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det0_nt2 = "df541a66a06983f90c60bc21fe4cc317",
    reg_res_trees_det1_nt1 = "ef724b6ac5ab79d23d13b6031fda14e0",
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
    pur_reg_components = "d839561fa1f77e429ab185db3711ef50",
    pur_reg_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_predict = "535a2cd73cf14947071df34e20eaab11",
    pur_reg_nt2_components = "d839561fa1f77e429ab185db3711ef50",
    pur_reg_nt2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_nt2_predict = "535a2cd73cf14947071df34e20eaab11",
    pur_reg_mode1_components = "30432c23c98b6558ef55e04879542d36",
    pur_reg_mode1_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_mode1_predict = "481c2302637dffaa23e872947fb6fedc",
    pur_reg_cap2_components = "db5778d6ea4116709de3be3ad5b518b9",
    pur_reg_cap2_names = "42622b48bfaa4d8b2d366e09a55c3899",
    pur_reg_cap2_predict = "f5847f7cf489166d19267149a735c6f1",
    pur_cls_components = "7d94825db5b12d04e8c96a52f7979d0f",
    pur_cls_names = "794adb1c06fff2ed25b9afca6d9c654c",
    pur_cls_predict = "1ec2de9b14b57d2abcda9fc799c525dd",
    pur_cls3_components = "3caa5837eac094639ac104e99751fdd8",
    pur_cls3_names = "28a2f21ce2d5f51e5a2bb15a2eee3069",
    pur_cls3_predict = "ddd5f951fc06622f6dcf0c4ecf461bfb"
  )
  got <- port_fit_digests()
  expect_identical(names(got), names(expected))
  for (nm in names(expected)) {
    expect_identical(got[[nm]], expected[[nm]], label = nm)
  }
})
