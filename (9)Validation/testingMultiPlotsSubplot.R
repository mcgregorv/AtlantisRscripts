##  Split the screen into two rows and one column, defining screens 1 and 2 and 3.
split.screen( figs = c( 3, 1 ) )

##  Split screen 1 into one row and three columns, defining screens 4, and 5.
split.screen( figs = c( 1, 2 ), screen = 1 )

##  Split screen 2 into one row and two columns, defining screens 6 and 7.
split.screen( figs = c( 1, 2 ), screen = 2 )

##  Split screen 3 into one row and two columns, defining screens 8 and 9.
split.screen( figs = c( 1, 2 ), screen = 3 )


##  The first plot is located in screen 3:
screen( 4 )
plot( rnorm( n = 10 ), col = "red", main = "plot 1" )
subplot(plot(rep(1,2), type="l"), x=4,y=1)

screen( 5 )
plot( rnorm( n = 10 ), col = "red", main = "plot 1" )



######################################

##  Split the screen into two rows and one column, defining screens 1 and 2.

split.screen( figs = c( 2, 1 ) )

##  Split screen 1 into one row and three columns, defining screens 3, 4, and 5.

split.screen( figs = c( 1, 3 ), screen = 1 )

##  Split screen 2 into one row and two columns, defining screens 6 and 7.

split.screen( figs = c( 1, 2 ), screen = 2 )

##  The first plot is located in screen 3:

screen( 3 )
plot( rnorm( n = 10 ), col = "red", main = "plot 1" )
subplot(plot(rep(1,2), type="l"), x=4,y=1)

##  The second plot is located in screen 4:

screen( 4 )
plot( runif( n = 10 ), col = "blue", main = "plot 2" )

##  The third plot is located in screen 5:

screen( 5 )
plot( rt( n = 10, df = 8 ), col = "springgreen4", main = "plot 3" )

##  The fourth plot is located in screen 6:

screen( 6 )
plot( rpois( n = 10, lambda = 2 ), col = "black", main = "plot 4" )

##  The fifth plot is located in screen 7:

screen( 7 )
plot( rf( n = 10, df1 = 4, df2 = 8 ), col = "gray30", main = "plot 5" )

##  Close all screens.

close.screen( all = TRUE )
