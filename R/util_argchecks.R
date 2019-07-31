.check_type_in_ami_nuc_unt <- function(type) {
  if (!type %in% c("unt", "ami", "nuc")) 
    stop("'type' has to be one of 'unt', 'ami', 'nuc'")
}

.check_nc_type_in_ami_nuc <- function(type) {
  if (missing(type) ||
      !type %in% c("nuc", "ami"))
    stop("in no_check mode 'type' needs to be one of 'nuc', 'ami'")
}

.check_sqstr_proper_char <- function(sq) {
  if (!is.character(sq) ||
      length(sq) == 0 ||
      any(is.na(sq)) ||
      any(is.null(sq))) 
    stop("'sq' has to be a positive-length vector of non-NA, non-NULL strings")
}

.check_is_clean_in_TRUE_FALSE_NULL <- function(is_clean) {
  if (!is.null(is_clean) &&
      !is_clean %in% c(TRUE, FALSE)) 
    stop("'is_clean' has to be TRUE, FALSE or NULL")
}

.check_nc_is_clean_in_TRUE_FALSE <- function(is_clean) {
  if (!is_clean %in% c(TRUE, FALSE)) 
    stop("'is_clean' has to be TRUE, FALSE in no_ceck mode")
}

.check_inds_are_numeric <- function(indices) {
  if (!(is.numeric(indices) && 
        floor(indices) == indices)) {
    stop("'indices' has to be an integer vector")
  }
}

.check_file_is_char <- function(file) {
  if (!is.character(file) ||
      !(length(file) == 1)) {
    stop("'file' has to be a string giving file to read from")
  }
}