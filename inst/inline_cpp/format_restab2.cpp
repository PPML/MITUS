// Using the inline package, we define the body of a function in a .cpp file here 
// and specify its arguments when calling it with cxxfunction.

// Convert ResTab, ResTabus, ResTabfb to NumericVectors 
NumericVector rt(ResTab);
NumericVector rtus(ResTabus);
NumericVector rtfb(ResTabfb);

// Record the dimensions of ResTab, ResTabus, ResTabfb
NumericVector dim = rt.attr("dim");
NumericVector dim2 = rtus.attr("dim");
NumericVector dim3 = rtfb.attr("dim");

// Make a NumericMatrix version of res_tab2
// When res_tab2 comes into the C++ code, it should already be a 
// NumericMatrix, but always better to be explicit than implicit. 
NumericMatrix rt2(res_tab2);

// This try ensures that the dimensions of rt, rtus, rtfb agree.
// If not, raise an error.
try {
	// make sure order of arrays matches
	if (dim.size() != dim2.size() | dim2.size() != dim3.size()) {
		throw std::invalid_argument("Dimensions of ResTab, ResTabus, and ResTabfb must match.");
	}
	// make sure dimensions match
	for (int j = 0; j < dim.size(); j++) {
		if (dim[j] != dim2[j] | dim2[j] != dim3[j]) {
			throw std::invalid_argument("Dimensions of ResTab, ResTabus, and ResTabfb must match.");
		}
	}
} catch(std::exception &ex) {	
	forward_exception_to_r(ex);
} catch(...) { 
	::Rf_error("c++ exception (unknown reason)"); 
}

// Construct sizes of increasing sequential dimensional subspaces,
// i.e. (unit (1), row-length, face-size, cube-size, etc...)
NumericVector dim_mults(dim.size());
for (int j = 0; j < dim.size(); j++) {
	dim_mults[j] = 1;
	if (j != 0) {
		for (int k = 0; k < j; k++) {
			dim_mults[j] = dim_mults[j] * dim[k];
		}
	}
}


// Initialize variables with the proper size to contain: 
// - ii: 
//      a row of the rt table
//
// - i: 
//      a version of ii, not containing the empty last cell which we'll be
//      computing
//
// - rt_idx: 
//     the position (in 4 dimensions) in ResTab, ResTabus, ResTabfb of the data
//     that we're using to fill in rt2
//
// - idx:
//     the position (in 1 dimension) in rt, rtus, or rtfb of the data that
//     we're going to use to fill in rt2
//
NumericVector ii(7);
IntegerVector i(8);
IntegerVector rt_idx(4);
int idx;

// Main logic loop:
//
// - Iterate through the rows of rt2 (which has yet to be filled in with real data) 
//
// - Compute from the row's entries how to look up the corresponding data from
//   rt rtus, or rtfb accordingly
//
// - Enter the new value into the last column of rt2
//
for (int j = 0; j < rt2.nrow(); j++) {
	// Get the jth row of rt2
	// The columns are: 
	//   0 - outcome:
	//         ltbi_000s, 
	//         pct_ltbi, 
	//         tb_incidence_000s, 
	//         tb_incidence_per_mil, 
	//         tb_mortality_000s, 
	//         tb_mortality_per_mil
	//   1 - scenario:
	//         base_case
	//   2 - population:
	//         all_populations,
	//         usb_population,
	//         fb_population
	//   3 - (small) age groups:
	//         0-4, 5-14, 15-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75-84, 85-94, 95+
	//   4 - comparator:
	//         absolute_value,
	//         pct_basecase_same_year,
	//         pct_basecase_2018
	//   5 - year:
	//         2018-2049
	//   6 - statistic:
	//         mean, ci_high_ci_low
	//   7 - value: 
	//         empty, filled in by the following loop
	//
	ii = rt2(j, _);
	// Turn ii into an IntegerVector for comparison with integers
	i = IntegerVector::create(ii[0],ii[1],ii[2],ii[3],ii[4],ii[5],ii[6]);


	// Outermost if-else Switches Arithmetic Based on Comparator in i[4]
	
	if (i[4] == 1) { 
	// ################## Report values in Absolute value ##################

		if (i[2] == 1) { 
			// -------------------- All Population -- Absolute Values --------------------

			// Determine 4-dimensional data position in ResTab
			IntegerVector rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1); 
			idx = 0; // Construct 1-dimensional data position in rt
			for (int k = 0; k < dim.size(); k++) { 
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rt[idx]; // Fill in data from rt

		} else if (i[2] == 2) { 
			// -------------------- US-Born Population -- Absolute Values --------------------

			// Determine 4-dimensional data position in ResTabus
			IntegerVector rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtus
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtus[idx]; // Fill in data from rtus


		} else if (i[2] == 3) { 
			// -------------------- Non-US-Born Population -- Absolute Values --------------------

			// Determine 4-dimensional data position in ResTabfb
			IntegerVector rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtfb
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtfb[idx]; // Fill in data from rtfb


		} else { 
			// -------------------- Throw an error if i[2] is not 1, 2, or 3 --------------------
			Rcpp::Rcout << "The population " << i[2] << " was referenced in res_tab2.";
			throw std::range_error("The population referenced in res_tab2 has not been defined in C++");
		}


	} else if (i[4] == 2) { 
		// ################## Report Values as Percent of the Basecase in the Same Year ###################

		if (i[2] == 1) { // All
			// ------------------ All Population -- Percent Basecase in Same Year ------------------

			// Determine 4-dimensional numerator data position in ResTab
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rt
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rt[idx]; // Fill in numerator data from rt

			// Determine 4-dimensional denominator data position in ResTab
			rt_idx = IntegerVector::create(0,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rt
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			} 

			// Divide by denominator data from rt and format as percent
			rt2(j,6) = rt2(j,6) / rt[idx] * 100; 


		} else if (i[2] == 2) { 
			// ------------------ US-Born Population -- Percent Basecase in Same Year -----------------

			// Determine 4-dimensional numerator data position in ResTabus
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtus
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtus[idx]; // Fill in numerator data from rtus

			// Determine 4-dimensional denominator data position in ResTabus
			rt_idx = IntegerVector::create(0,ii[5]-1,ii[0],ii[3]-1);
			idx = 0;
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}

			// Divide by denominator data from rtus and format as percent
			rt2(j,6) = rt2(j,6) / rtus[idx] * 100;

		} else if (i[2] == 3) { 
			// ------------------ Non-US-Born Population -- Percent Basecase in Same Year ------------------

			// Determine 4-dimensional numerator data position in ResTabfb
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtfb
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtfb[idx]; // Fill in numerator data from rtfb

			// Determine 4-dimensional denominator data position in ResTabfb
			rt_idx = IntegerVector::create(0,ii[5]-1,ii[0],ii[3]-1);
			idx = 0;
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}

			// Divide by denominator data from rtfb and format as percent
			rt2(j,6) = rt2(j,6) / rtfb[idx] * 100;

		} else {
			// ------------------------ Throw an error if i[2] is not 1, 2, or 3 ------------------------
			Rcpp::Rcout << "The population " << i[2] << " was referenced in res_tab2.";
			throw std::range_error("The population referenced in res_tab2 has not been defined in C++");
		}


	} else if (i[4] == 3) { 
		// ################## Report Values as Percent of the Basecase in 2018 ###################
		
		if (i[2] == 1) { 
			// ------------------ All Population -- Percent Basecase in 2018 ------------------------
			
			// Determine 4-dimensional numerator data position in ResTab
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rt
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rt[idx]; // Fill in numerator data from rt

			// Determine 4-dimensional denominator data position in ResTab
			rt_idx = IntegerVector::create(0,ii[5]-1,0,ii[3]-1);
			idx = 0;
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}

			// Divide by denominator data from rt and format as percent
			rt2(j,6) = rt2(j,6) / rt[idx] * 100;

		} else if (i[2] == 2) {
			// ------------------ US-Born Population -- Percent Basecase in 2018 ------------------------

			// Determine 4-dimensional numerator data position in ResTabus
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtus
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtus[idx]; // Fill in numerator data from rtus

			// Determine 4-dimensional denominator data position in ResTabus
			rt_idx = IntegerVector::create(0,ii[5]-1,0,ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtus
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}

			// Divide by denominator data from rtus and format as percent
			rt2(j,6) = rt2(j,6) / rtus[idx] * 100;


		} else if (i[2] == 3) { 
			// ------------------ Non-US-Born Population -- Percent Basecase in 2018 ------------------------
			
			// Determine 4-dimensional numerator data position in ResTabfb
			rt_idx = IntegerVector::create(ii[1]-1,ii[5]-1,ii[0],ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtfb
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}
			rt2(j,6) = rtfb[idx]; // Fill in numerator data from rtfb

			// Determine 4-dimensional denominator data position in ResTabfb
			rt_idx = IntegerVector::create(0,ii[5]-1,0,ii[3]-1);
			idx = 0; // Construct 1-dimensional data position in rtfb
			for (int k = 0; k < dim.size(); k++) {
				idx = idx + dim_mults[k] * rt_idx[k];
			}

			// Divide by denominator data from rtus and format as percent
			rt2(j,6) = rt2(j,6) / rtfb[idx] * 100;

		} else {
			// ------------------ Throw an error if i[2] is not 1, 2, or 3 ------------------------
			Rcpp::Rcout << "The population " << i[2] << " was referenced in res_tab2.";
			throw std::range_error("The population referenced in res_tab2 has not been defined in C++");
		}
	} else {
		// ------------------ Throw an error if i[3] is not 1, 2, or 3 ------------------------
		Rcpp::Rcout << "The comparator " << i[4] << " was referenced in res_tab2.";
		throw std::range_error("The comparator referenced in res_tab2 has not been defined in C++");
	}
}

return rt2;
