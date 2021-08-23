OBJECTS1 = delete_CR2 avgPlain fractal_dim2 fractal_dim_mass complexityBrinkhoff complexity1 complexity3 complexity5 complexityPage convexification convexity convexity_jovisa2 k_regularity match_poly_hausdorf normalise_shape_mom3 pix_transf ramer super2pixel avg12

all: $(OBJECTS1)

clean:
	rm -f $(OBJECTS1)

$(OBJECTS1): $(OBJECTS1:=.c)
	cc -o $@ $@.c -lm

