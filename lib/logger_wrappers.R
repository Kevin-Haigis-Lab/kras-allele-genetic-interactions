
# Wrapper functions for common logging patterns

# Check that there are rows in a data frame `df`
log_rows <- function(logger_obj, df, df_name) {
    if (nrow(df) > 0) {
        info(logger_obj, paste0("`", df_name, "` has ", nrow(df), " rows."))
    } else {
        error(logger_obj, paste0("`", df_name, "` has 0 rows."))
    }
}
