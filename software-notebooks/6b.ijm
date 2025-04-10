#@ File (label = "Input directory", style = "directory") par_input
#@ File (label = "Output directory", style = "directory") par_output
#@ String (label = "File suffix", value = ".czi") par_suffix
#@ Double (label = "Surface Threshold", value = 15000) par_surface_threshold
#@ Double (label = "Lipid Threshold", value = 1500) par_lipid_threshold
#@ Integer (label = "Phase slice in Focus", value = 1) par_phase_focus_slice
#@ Integer (label = "Lipid slice in Focus", value = 3) par_lipid_focus_slice

//----------------------------------------------------------------//
// Quantification of the area occupied by cells in a phase image  //
// Quantification of the area occupied by lipid droplets          //
//                                                                //
// (c) 2023 R.Harkes & B.vd.Broek - NKI                           //
// Contact: r.harkes@nki.nl                                       //
// Version 1.0                                                    //
//----------------------------------------------------------------//

table_name = "Measurement_Table";  
Table.create(table_name);
Table.reset(table_name);
run("Clear Results");
run("Conversions...", "scale");
setBatchMode(true);
processFolder(par_input, par_output, par_suffix);
setBatchMode(false);
// Save table
tablePath = par_output + File.separator + "Result.csv";
print(tablePath);
Table.save(tablePath, table_name);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input, output, suffix) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		file = input  + File.separator + list[i];
		if(File.isDirectory(file)){
			outFolder = output + File.separator + list[i];
			File.makeDirectory(outFolder);
			processFolder(file, outFolder, suffix);
		}
		if(endsWith(list[i], suffix)){
			processFile(file, output, suffix);
		}
	}
}

function processFile(fileIn, outFolder, suffix) {
	print("Processing: " + fileIn);
	filename = File.getName(fileIn);
	filename_no_suffix = substring(filename, 0, filename.length - suffix.length);
	run("Close All");
	run("Bio-Formats Importer", "open=[" + fileIn + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	rename(filename_no_suffix);
	
	// Phase image analysis
	Stack.setPosition(3, par_phase_focus_slice, 1);
	run("Duplicate...", " ");
	phase_slice = "phaseslice";
	rename(phase_slice);
	run("Duplicate...", " ");
	run("Gaussian Blur...", "sigma=75");
	rename(phase_slice+"_blur");
	imageCalculator("Subtract create 32-bit", phase_slice, phase_slice+"_blur");
	rename(phase_slice+"_corrected");
	run("Enhance Contrast", "saturated=0.35");
	run("Duplicate...", " ");
	run("Variance...", "radius=20");
	setThreshold(par_surface_threshold, 1e30);
	cell_area = getValue("Area limit");
	curr_row = Table.size(table_name);
	Table.set("File Dir", curr_row, File.getDirectory(fileIn), table_name);
	Table.set("File Name", curr_row, filename_no_suffix, table_name);
	Table.set("Cell Area (um2)", curr_row, cell_area, table_name);
	run("Convert to Mask");
	rename(phase_slice+"_mask");
	
	// Lipid droplet detection
	selectWindow(filename_no_suffix);
	Stack.setPosition(2, par_lipid_focus_slice, 1);
	run("Duplicate...", " ");
	drop_slice = "dropslice";
	rename(drop_slice);
	setThreshold(par_lipid_threshold, 65535);
	green_area = getValue("Area limit");
	Table.set("Lipid Area (um2)", curr_row, green_area, table_name);
	Table.set("% Area (um2)", curr_row, green_area/cell_area*100, table_name);
	Table.update;
	selectWindow(drop_slice);
	run("Convert to Mask");
	rename(drop_slice+"_mask");
	
	// Save corrected phase image with background overlay in magenta and lipid droplets in green
	selectWindow(phase_slice+"_mask");
	run("Invert");
	run("Magenta");
	selectWindow(phase_slice+"_corrected");
	run("Add Image...", "image=" +phase_slice+"_mask x=0 y=0 opacity=50 zero");
	selectWindow(drop_slice+"_mask");
	run("Green");
	selectWindow(phase_slice+"_corrected");
	run("Add Image...", "image=" +drop_slice+"_mask x=0 y=0 opacity=50 zero");
	run("Flatten");
	fileOut = outFolder + File.separator + filename_no_suffix + "_phase.jpg";
	saveAs("jpg", fileOut);
}