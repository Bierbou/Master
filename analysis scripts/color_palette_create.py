import enve
nlevels = 21
infile = r"C:\Users\kevin\Documents\cropped_image.png"
outfile = "example_palette.cpal"

img = enve.image()
err=img.load(infile)
size = img.dims
x = 0
total_y = int(size[1])
incr = total_y / (nlevels - 1)
rgb_val = []

for y in range(0,total_y - incr,incr):
	rgb = img.getpixel(x,y)
	rval = '%5.4f' % (float(rgb[0]) / 255.0)
	gval = '%5.4f' % (float(rgb[1]) / 255.0)
	bval = '%5.4f' % (float(rgb[2]) / 255.0)
	if y % incr == 0:
		rgb_val.append([rval,gval,bval])
# Now, little fiddly to put last value in the rgb_val too:
rgb = img.getpixel(x,total_y - 1)
rval = '%5.4f' % (float(rgb[0]) / 255.0)
gval = '%5.4f' % (float(rgb[1]) / 255.0)
bval = '%5.4f' % (float(rgb[2]) / 255.0)
rgb_val.append([rval,gval,bval])

# Write out the cpal file
ofile = open(outfile,"w")
ofile.write("number_of_levels   " + str(len(rgb_val)) + "\n")
ofile.write("colors\n")
for i in range(len(rgb_val)):
	ofile.write(str(rgb_val[i][0]) + " " + str(rgb_val[i][1]) + " " + str(rgb_val[i][2]) + "\n")
ofile.close()

