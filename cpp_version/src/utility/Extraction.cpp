#define __USE_MINGW_ANSI_STDIO 0
#include"cnpy.h"
#include"Extraction.h"


int get_offset_regions(std::string offset_regions_filepath, std::unique_ptr<int>& offset_region_ptr) {
	// assign the ptr to the data to offset_region_ptr
	// return the length the data ( the second dimension)
	cnpy::NpyArray offset_region_array = cnpy::npy_load(offset_regions_filepath);
	assert(offset_region_array.shape.size() == 2);
	assert(offset_region_array.shape[0] == 8);
	offset_region_ptr = std::make_unique<int>(offset_region_array.data<int>());
	return offset_region_array.shape[1];
}

void make_features(std::string image_filepath, std::unique_ptr<int> offset_region_ptr, int num_comparisons)
{
	//load it into a new array

	//const size_t dim[3] = { IMAGE_HEIGHT, IMAGE_WIDTH, IMAGE_CHANNELS };

	cnpy::NpyArray image_integral_array = cnpy::npy_load(image_filepath);

	assert(image_integral_array.shape.size() == 3);
	/*for (size_t i = 0; i < image_integral_array.shape.size(); ++i) {
		assert(dim[i] == image_integral_array.shape[i]);
	}*/

	double* loaded_data = image_integral_array.data<double>();

	const size_t height = image_integral_array.shape[0];
	const size_t width = image_integral_array.shape[1];
	const size_t num_channels = image_integral_array.shape[2];

	//for (size_t i = 0; i < height; ++i) {
	//	for (size_t j = 0; j < width; ++j) {
	//		int	x1 = i + int(u1[m] / depth)
	//		int	y1 = j + int(v1[m] / depth)
			//	x2 = i + int(u2[m] / depth)
			//	y2 = j + int(v2[m] / depth)
			//	w1n = int(w1[m] / depth)
			//	h1n = int(h1[m] / depth)
			//	w2n = int(w2[m] / depth)
			//	h2n = int(h2[m] / depth)
			//	R1 = float(w1n * h1n)
			//	R2 = float(w2n * h2n)
	//	}
	//}


	//for m in range(n_comparisons) :
	//	x1 = i + int(u1[m] / depth)
	//	y1 = j + int(v1[m] / depth)
	//	x2 = i + int(u2[m] / depth)
	//	y2 = j + int(v2[m] / depth)
	//	w1n = int(w1[m] / depth)
	//	h1n = int(h1[m] / depth)
	//	w2n = int(w2[m] / depth)
	//	h2n = int(h2[m] / depth)
	//	R1 = float(w1n * h1n)
	//	R2 = float(w2n * h2n)
	//	# design feature
	//	# offset - difference as feature
	//	data = numpy.zeros(n_channels, dtype = float)
	//	for k in range(n_channels) :
	//		try :
	//		data[k] = (float(summed_area[x1 + w1n, y1 + h1n, k] + summed_area[x1, y1, k] \
	//			- summed_area[x1, y1 + h1n, k] - summed_area[x1 + w1n, y1, k]) / R1) \
	//		- (float(summed_area[x2 + w2n, y2 + h2n, k] + summed_area[x2, y2, k]
	//			- summed_area[x2, y2 + h2n, k] - summed_area[x2 + w2n, y2, k]) / R2)
	//		except IndexError :
	//data[k] = LARGE_VAL
	//	data_entry.extend(data)



	return;
}
