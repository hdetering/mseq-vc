require(UpSetR)

plot_upset <- function( df, lbl_sets ) {
  upset(df, sets = lbl_sets, 
        keep.order = TRUE, 
        order.by = 'freq',
        mb.ratio = c(0.5, 0.5),
        number.angles = 45)
}