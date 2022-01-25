library(GeoplotR)

pdf('test.pdf')

# test 1:
xytest('Template.json',polygons=FALSE,short=TRUE)
xytest('Template.json',polygons=FALSE,short=FALSE)
xytest('Template.json',polygons=TRUE,short=TRUE)

# test 2:
xyztest('TiZrY.json',polygons=TRUE,short=TRUE)
xyztest('TiZrY.json',polygons=TRUE,short=FALSE)
xyztest('TiZrY.json',polygons=FALSE,short=TRUE)

# test 3:
xytest('Cr_Y.json',polygons=FALSE)
xytest('Cr_Y.json',log='x',polygons=FALSE)
xytest('Cr_Y.json',log='xy',polygons=FALSE)
xytest('Cr_Y.json',log='xy',polygons=FALSE,smooth=TRUE)
xytest('Cr_Y.json',log='xy',polygons=TRUE)

dev.off()
