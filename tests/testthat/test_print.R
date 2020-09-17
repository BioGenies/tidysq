# # SETUP ----
# sq_dna <- construct_sq_dna(c("ACTGTC", "CGCGTTA"), is_clean = TRUE)
# sq_rna <- random_sq(25, 10, "rna", is_clean = FALSE)
# sq_ami <- construct_sq_ami(c("APOPNIQEV", "CSVMIBF"), is_clean = FALSE)
# sq_unt <- construct_sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"))
# sq_atp <- construct_sq(c("AmiDnaDnaAmi", "RnaDnaAmiDna"),
#                        non_standard = c("Ami", "Dna", "Rna"))
# sq_enc <- encode(sq_dna, c(A = 1, C = 2, G = 3, T = 21))
# sq_empty <- construct_sq_rna(character(), is_clean = TRUE)
# 
# header_regexp <- "\\w{3} \\(.+?\\)(, cln \\(cleaned\\))? sequences list[:|( of length 0)]"
# footer_regexp <- function(sq_length, option_length) {
#   if (sq_length <= option_length) {
#     # NA value is used by expect_output to signify no output
#     NA
#   } else {
#     sprintf("printed %s out of %s", option_length, sq_length)
#   }
# }
# 
# # HEADER PRINTING ----
# test_that("obj_print_header() satisfies header regexp", {
#   print("\n")
#   obj_print_header(sq_rna)
#   print("\n")
#   print(sq_rna)
#   print("\n")
#   expect_output(capture_output(obj_print_header(sq_dna)), regexp = header_regexp)
#   expect_output(obj_print_header(sq_rna), regexp = header_regexp)
#   expect_output(obj_print_header(sq_ami), regexp = header_regexp)
#   expect_output(obj_print_header(sq_unt), regexp = header_regexp)
#   expect_output(obj_print_header(sq_atp), regexp = header_regexp)
#   expect_output(obj_print_header(sq_enc), regexp = header_regexp)
#   expect_output(obj_print_header(sq_empty), regexp = header_regexp)
# })
# 
# # FOOTER PRINTING ----
# test_that("obj_print_footer() satisfies footer regexp for default number of sequences", {
#   # Tests with default value of 10 sequences
#   withr::local_options(list(tidysq_p_max_sequences = 10))
#   expect_output(obj_print_footer(sq_dna),
#                 regexp = footer_regexp(length(sq_dna), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_rna),
#                 regexp = footer_regexp(length(sq_rna), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_ami),
#                 regexp = footer_regexp(length(sq_ami), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_unt),
#                 regexp = footer_regexp(length(sq_unt), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_atp),
#                 regexp = footer_regexp(length(sq_atp), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_enc),
#                 regexp = footer_regexp(length(sq_enc), getOption("tidysq_p_max_sequences")))
#   expect_output(obj_print_footer(sq_empty),
#                 regexp = footer_regexp(length(sq_empty), getOption("tidysq_p_max_sequences")))
# })
