#!/bin/csh -f

# store temporary outputs here
if !(-d tmp) then
	mkdir tmp
endif

# delete any old results
foreach i (fractal_dim2.dat fractal_dim_massRobust.dat complexity3NormG16.dat complexityBrinkhoff.dat complexity1.dat complexity5NormG16.dat complexityPage.dat convexA.dat convexJ.dat convexL.dat convexificationArea.dat convexificationRevArea.dat k_regularity.dat)
	if -e $i then
		rm $i
	endif
end

# compute measures for each view
foreach i (*.pix)
	### some complexity measures

	./complexityBrinkhoff $i >> complexityBrinkhoff.dat
	./complexity1 $i 4 | grep -v reset >> complexity1.dat
	./complexityPage $i >> complexityPage.dat

	./normalise_shape_mom3 -i $i -o ./tmp/norm -a 100000 -g 1 > /dev/null
	./complexity3 ./tmp/norm 16 >> complexity3NormG16.dat
	./complexity5 ./tmp/norm 16 >> complexity5NormG16.dat

	### a "regularity" measure

	./k_regularity $i >> k_regularity.dat

	### some convexity measures

	./ramer -i $i -o ./tmp/lines -d 0.5
	./super2pixel -i ./tmp/lines -o ./tmp/pixels

	./convexification -i ./tmp/pixels -o ./tmp/convexification > ./tmp/convexification2
	grep area ./tmp/convexification2 >> convexificationArea.dat

	./convexification -r -i ./tmp/pixels -o ./tmp/convexification > ./tmp/convexification2
	grep area ./tmp/convexification2 >> convexificationRevArea.dat

	./convexity $i > ./tmp/ccc
	grep area ./tmp/ccc >> convexA.dat
	grep length ./tmp/ccc >> convexL.dat

	# use the standard level of polygonalisation that I've used in the past for this measure
	./ramer -i $i -o ./tmp/tmp1 > /dev/null
	./convexity_jovisa2 -i ./tmp/tmp1 | grep "convexity =" >> convexJ.dat

	### some fractal measures

	./pix_transf -i $i -o ./tmp/pixF -F
	./fractal_dim2 ./tmp/pixF | grep dimension >> fractal_dim2.dat
	./fractal_dim_mass -i ./tmp/pixF -o ./tmp/fractal -I 1 -r ; grep dimension ./tmp/fractal >> fractal_dim_massRobust.dat
end

# strip leading text from .dat files
./run_keep_end

# average scores across 12 views and display results
foreach i (fractal_dim2.dat fractal_dim_massRobust.dat complexity3NormG16.dat complexityBrinkhoff.dat complexity1.dat complexity5NormG16.dat complexityPage.dat convexA.dat convexJ.dat convexL.dat convexificationArea.dat convexificationRevArea.dat k_regularity.dat)
	echo === $i
	./avg12 < $i
end

### method based on similarity across views
# separately run this method since it involves performing a pairwise comparison of all the views for a given model
# unlike the other methods that operate on each view independently
# it also means that this method is slow!
./run_complexityView1.sh
# display results
echo === complexityView1.dat
cat complexityView1.dat
