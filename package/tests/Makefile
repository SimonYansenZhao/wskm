update:
	R CMD BATCH ewkm01.R
	mv ewkm01.Rout ewkm01.Rout.save
	rm .RData
	R CMD BATCH fgkm_test.R
	mv fgkm_test.Rout fgkm_test.Rout.save
	R CMD BATCH twkm_test.R
	mv twkm_test.Rout twkm_test.Rout.save
