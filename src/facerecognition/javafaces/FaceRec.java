/**
 * Changes inspired from
 * http://jjandroidlearning.blogspot.ro/2012/04/change-from-javafaces-into-android.html
 */

package facerecognition.javafaces;

import java.io.*;
import java.util.*;
import java.util.logging.Handler;

import android.annotation.SuppressLint;
import android.content.Context;
import android.graphics.Bitmap;
import android.graphics.BitmapFactory;
import android.graphics.Canvas;
import android.graphics.Color;
import android.graphics.ColorMatrix;
import android.graphics.ColorMatrixColorFilter;
import android.graphics.Paint;
import android.util.Log;

import facerecognition.utils.Utils;
import facerecognition.utils.ValueIndexPair;

@SuppressLint("DefaultLocale")
public class FaceRec {
	private Context androidContext = null;
	private final static String LOG_TAG = "FaceRec";
	
	private FaceBundle bundle;
	private double[][] weights;
	Handler fh;

	public FaceRec(Context ctx) {
		this.androidContext = ctx;
	}

	public MatchResult findMatchResult(Bitmap image, int selectedeigenfaces,
			double thresholdVal) {

		boolean match = false;
		String message = null;
		String matchingFileName = "";
		double minimumDistance = 0.0;
		try {
			Log.d(LOG_TAG, "findMatchResult+");
			
			checkImageSizeCompatibility(image);
			
			Matrix2D inputFace = getNormalisedInputFace(image);
			Log.d(LOG_TAG, "if:" + inputFace.get(0, 0) + " " + inputFace.get(0, 1));
			
			inputFace.subtract(new Matrix2D(bundle.getAvgFace(), 1));
			Log.d(LOG_TAG, "if:" + inputFace.get(0, 0) + " " + inputFace.get(0, 1));

			Matrix2D inputWts = getInputWeights(selectedeigenfaces, inputFace);
			Log.d(LOG_TAG, "iw rows:" + inputWts.rows());
			Log.d(LOG_TAG, "iw cols:" + inputWts.columns());
			Log.d(LOG_TAG, "iw:" + inputWts.get(0, 0));

			double[] distances = getDistances(inputWts);
			Log.d(LOG_TAG, "distances:");
			for (int i = 0; i < distances.length; i++) {
				Log.d(LOG_TAG, ""+distances[i]);
			}
			
			ImageDistanceInfo distanceInfo = getMinimumDistanceInfo(distances);
			minimumDistance = Math.sqrt(distanceInfo.getValue());
			matchingFileName = getMatchingFileName(distanceInfo);

			if (minimumDistance > thresholdVal) {
				message = "no match found, try higher threshold";
				Log.d("FaceRec", message);
			} else {
				match = true;
				message = "matching image found";
				Log.d("FaceRec", message);
				Log.d("FaceRec", minimumDistance + " vs " + thresholdVal);

			}
		} catch (Exception e) {
			Log.e("FaceRec", "Exception thrown..");
			e.printStackTrace();
			return new MatchResult(false, "", Double.NaN, e.getMessage());
		}
		
		Log.d(LOG_TAG, "findMatchResult- OK");
		return new MatchResult(match, matchingFileName, minimumDistance, message);
	}

	public MatchResult findMatchResult(String imageFileName,
			int selectedeigenfaces, double thresholdVal) {
		boolean match = false;
		String message = null;
		String matchingFileName = "";
		double minimumDistance = 0.0;
		try {
			Log.d(LOG_TAG, "findMatchResult+ (String)");

			checkImageSizeCompatibility(imageFileName);
			
			Matrix2D inputFace = getNormalisedInputFace(imageFileName);
			Log.d(LOG_TAG, "if:" + inputFace.get(0, 0) + " " + inputFace.get(0, 1));
			
			inputFace.subtract(new Matrix2D(bundle.getAvgFace(), 1));
			Log.d(LOG_TAG, "if:" + inputFace.get(0, 0) + " " + inputFace.get(0, 1));

			Matrix2D inputWts = getInputWeights(selectedeigenfaces, inputFace);
			Log.d(LOG_TAG, "iw rows:" + inputWts.rows());
			Log.d(LOG_TAG, "iw cols:" + inputWts.columns());
			Log.d(LOG_TAG, "iw:" + inputWts.get(0, 0));

			double[] distances = getDistances(inputWts);
			Log.d(LOG_TAG, "distances:");
			for (int i = 0; i < distances.length; i++) {
				Log.d(LOG_TAG, ""+distances[i]);
			}
			
			ImageDistanceInfo distanceInfo = getMinimumDistanceInfo(distances);
			minimumDistance = Math.sqrt(distanceInfo.getValue());
			matchingFileName = getMatchingFileName(distanceInfo);

			if (minimumDistance > thresholdVal) {
				message = "no match found, try higher threshold";
			} else {
				match = true;
				message = "matching image found";
			}
		} catch (Exception e) {
			Log.d(LOG_TAG, "findMatchResult-e");
			return new MatchResult(false, "", Double.NaN, e.getMessage());
		}
		
		Log.d(LOG_TAG, "findMatchResult-");
		return new MatchResult(match, matchingFileName, minimumDistance, message);
	}

	private void checkImageSizeCompatibility(Bitmap image) throws IOException,
			FaceRecError {
		int height = image.getHeight();
		int width = image.getWidth();

		int facebundleWidth = this.bundle.getImageWidth();
		int facebundleHeight = this.bundle.getImageHeight();

		if ((height != facebundleHeight) || (width != facebundleWidth)) {
			throw new FaceRecError(
					"selected image dimensions does not match dimensions of other images");
		}

	}

	private void checkImageSizeCompatibility(String fileName) throws IOException,
			FaceRecError {
		Bitmap image = BitmapFactory.decodeFile(fileName);
		int height = image.getHeight();
		int width = image.getWidth();

		int facebundleWidth = this.bundle.getImageWidth();
		int facebundleHeight = this.bundle.getImageHeight();

		if ((height != facebundleHeight) || (width != facebundleWidth)) {
			throw new FaceRecError(
					"selected image dimensions does not match dimensions of other images");
		}

	}

	private String getMatchingFileName(ImageDistanceInfo distanceInfo) {
		List<String> imageFileNames = this.bundle.getImageFileNamesList();
		String matchingFileName = imageFileNames.get(distanceInfo.getIndex());
		return matchingFileName;
	}

	private Matrix2D getInputWeights(int selectedeigenfaces, Matrix2D inputFace) {
		double[][] eigenFacesArray = this.bundle.getEigenFaces();
		Matrix2D eigenFacesMatrix = new Matrix2D(eigenFacesArray);
		Matrix2D eigenFacesMatrixPart = eigenFacesMatrix
				.getSubMatrix(selectedeigenfaces);
		Matrix2D eigenFacesMatrixPartTranspose = eigenFacesMatrixPart.transpose();
		Matrix2D inputWts = inputFace.multiply(eigenFacesMatrixPartTranspose);
		return inputWts;
	}

	private Matrix2D getNormalisedInputFace(String imageFileName)
			throws FaceRecError {
		double[] inputFaceData = getImageData(imageFileName);
		Matrix2D inputFace = new Matrix2D(inputFaceData, 1);
		// Log.i("FaceRec", "inputface");
		// Log.i("FaceRec", inputFace.toString());
		inputFace.normalise();
		// Log.i("FaceRec", "normalised inputface");
		// Log.i("FaceRec", inputFace.toString());
		return inputFace;
	}

	private Matrix2D getNormalisedInputFace(Bitmap image) throws FaceRecError {
		double[] inputFaceData = getImageData(image);

		Matrix2D inputFace = new Matrix2D(inputFaceData, 1);
		inputFace.normalise();

		return inputFace;
	}

	private ImageDistanceInfo getMinimumDistanceInfo(double[] distances) {
		double minimumDistance = Double.MAX_VALUE;
		int index = 0;
		for (int i = 0; i < distances.length; i++) {
			if (distances[i] < minimumDistance) {
				minimumDistance = distances[i];
				index = i;
			}
		}
		return new ImageDistanceInfo(distances[index], index);
	}

	// TODO: debug why temp[i][j] has a NaN value
	private double[] getDistances(Matrix2D inputWt) {
		Matrix2D tempWt = new Matrix2D(this.weights);
		Log.d("FaceRec", "weights: " + tempWt.toString());
		
		double[] inputWtData = inputWt.flatten();
		tempWt.subtractFromEachRow(inputWtData);
		tempWt.multiplyElementWise(tempWt);
		double[][] temp = tempWt.toArray();
		double[] distances = new double[temp.length];
		for (int i = 0; i < temp.length; i++) {
			double sum = 0.0;
			for (int j = 0; j < temp[0].length; j++) {
				sum += temp[i][j];
				Log.d("FaceRec", "temp[i][j]:" + temp[i][j]);
			}
			
			Log.d("FaceRec", "sum:" + sum);
			distances[i] = sum;
		}
		return distances;
	}

	// TODO: refactor such for bw pics only
	private double[] getImageData(String imageFileName) throws FaceRecError {
		double[] result = null;
		int[] pixels = null;

		// loading image
		Bitmap bmp = BitmapFactory.decodeFile(imageFileName);

		// getting pixels
		pixels = new int[bmp.getWidth() * bmp.getHeight()];
		bmp.getPixels(pixels, 0, bmp.getWidth(), 0, 0, bmp.getWidth(),
				bmp.getHeight());

		// converting to double
		result = new double[pixels.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (double) (Color.blue(pixels[i]));
		}

		return result;
	}

	private double[] getImageData(Bitmap bmp) throws FaceRecError {
		double[] result = null;
		int[] pixels = null;

		// getting pixels
		pixels = new int[bmp.getWidth() * bmp.getHeight()];
		bmp.getPixels(pixels, 0, bmp.getWidth(), 0, 0, bmp.getWidth(),
				bmp.getHeight());

		// converting to double
		result = new double[pixels.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = (double) (Color.blue(pixels[i]));
		}

		return result;
	}

	private void doCalculations(String dir, List<String> imglist,
			int selectedNumOfEigenFaces) throws FaceRecError, IOException {
		FaceBundle b = createFaceBundle(imglist);
		double[][] wts = calculateWeights(b, selectedNumOfEigenFaces);
		this.bundle = b;
		this.weights = wts;
		writeCache(dir, b);
	}

	private double[][] calculateWeights(FaceBundle b, int selectedNumOfEigenFaces) {
		// TODO: experimenting
		System.gc();
		
		Matrix2D eigenFaces = new Matrix2D(b.getEigenFaces());
		Matrix2D eigenFacesPart = eigenFaces.getSubMatrix(selectedNumOfEigenFaces);
		Matrix2D adjustedFaces = new Matrix2D(b.getAdjustedFaces());
		Matrix2D eigenFacesPartTr = eigenFacesPart.transpose();
		Matrix2D wts = adjustedFaces.multiply(eigenFacesPartTr);
		return wts.toArray();
	}

	public FaceBundle createFaceBundle(List<String> filenames)
			throws FaceRecError, IOException {
		Bitmap[] bufimgs = getGrayScaleImages(filenames);
		checkImageDimensions(filenames, bufimgs);

		Matrix2D imagesData = getNormalisedImagesData(bufimgs);

		double[] averageFace = imagesData.getAverageOfEachColumn();
		imagesData.adjustToZeroMean();
		// Log.d("FaceRec", "imagesData adjusted ToZeroMean");
		 Log.d("FaceRec", imagesData.columns() + "x" + imagesData.rows());
		// Log.d("FaceRec", "done imageData");

		EigenvalueDecomposition egdecomp = getEigenvalueDecomposition(imagesData);
		double[] eigenvalues = egdecomp.getEigenValues();
		double[][] eigvectors = egdecomp.getEigenVectors();

		 Log.d("FaceRec", "eigenvalues");
		 Log.d("FaceRec", new Matrix2D(eigenvalues, 1).toString());
		
		 Log.d("FaceRec", "eigvectors");
		 Log.d("FaceRec", new Matrix2D(eigvectors).toString());

		TreeSet<ValueIndexPair> pairList = getSortedPairs(eigenvalues, eigvectors);
		eigenvalues = getSortedVector(pairList);
		eigvectors = getSortedMatrix(eigvectors, pairList);

		// Log.d("FaceRec", "AFTER SORTING");
		// Log.d("FaceRec", "eigenvalues");
		// Log.d("FaceRec", new Matrix2D(eigenvalues, 1).toString());
		//
		// Log.d("FaceRec", "eigvectors");
		// Log.d("FaceRec", new Matrix2D(eigvectors).toString());

		Matrix2D eigenFaces = getNormalisedEigenFaces(imagesData, new Matrix2D(
				eigvectors));
		int imageWidth = bufimgs[0].getWidth();
		createEigenFaceImages(eigenFaces, imageWidth);
		int imageHeight = bufimgs[0].getHeight();
		FaceBundle b = new FaceBundle(filenames, imagesData.toArray(), averageFace,
				eigenFaces.toArray(), eigenvalues, imageWidth, imageHeight);
		return b;
	}

	public double[] getSortedVector(TreeSet<ValueIndexPair> pairSet) {
		double[] sortedVector = new double[pairSet.size()];
		ValueIndexPair[] viArray = pairSet.toArray(new ValueIndexPair[0]);
		for (int i = 0; i < pairSet.size(); i++) {
			sortedVector[i] = viArray[i].getVectorElement();
		}
		return sortedVector;
	}

	public double[][] getSortedMatrix(double[][] origmatrix,
			TreeSet<ValueIndexPair> pairSet) {
		int rows = pairSet.size();
		int cols = origmatrix[0].length;
		double[][] sortedMatrix = new double[rows][cols];
		ValueIndexPair[] viArray = pairSet.toArray(new ValueIndexPair[0]);
		// fill a 2D array using data from rows of original matrix
		for (int i = 0; i < pairSet.size(); i++) {
			sortedMatrix[i] = origmatrix[viArray[i].getMatrixRowIndex()];
		}

		return sortedMatrix;
	}

	public TreeSet<ValueIndexPair> getSortedPairs(double[] aVector,
			double[][] aMatrix) {
		TreeSet<ValueIndexPair> pairSet = createPairs(aVector, aMatrix);
		return pairSet;
	}

	public TreeSet<ValueIndexPair> createPairs(double[] aVector,
			double[][] aMatrix) {
		TreeSet<ValueIndexPair> pList = null;
		if (aVector.length != aMatrix.length) {
			printError("matrix rows don't match items in vector ");
		} else {
			pList = new TreeSet<ValueIndexPair>();
			for (int i = 0; i < aVector.length; i++) {
				ValueIndexPair dp = new ValueIndexPair(aVector[i], i);
				pList.add(dp);
			}
		}
		return pList;
	}

	private EigenvalueDecomposition getEigenvalueDecomposition(Matrix2D imagesData) {
		Matrix2D imagesDataTr = imagesData.transpose();
		Matrix2D covarianceMatrix = imagesData.multiply(imagesDataTr);
		EigenvalueDecomposition egdecomp = covarianceMatrix
				.getEigenvalueDecomposition();
		return egdecomp;
	}

	public void createEigenFaceImages(Matrix2D eigenfaces, int imgwidth)
			throws IOException {
		Log.i("FaceRec", "creating eigenfaces");
		double[][] eigenfacesArray = eigenfaces.toArray();
		String fldrname = this.getAbsoluteFolderPath("eigenfaces");
		makeNewFolder(fldrname);
		String prefix = "eigen";
		String ext = ".png";
		for (int i = 0; i < eigenfacesArray.length; i++) {
			double[] egface = eigenfacesArray[i];
			String filename = fldrname + File.separator + prefix + i + ext;
			createImageFromArray(filename, egface, imgwidth);
		}
		Log.i("FaceRec", "created eigenfaces.");
	}

	private Matrix2D getNormalisedEigenFaces(Matrix2D imagesData,
			Matrix2D eigenVectors) {
		Matrix2D eigenFaces = eigenVectors.multiply(imagesData);
		double[][] eigenFacesData = eigenFaces.toArray();
		for (int i = 0; i < eigenFacesData.length; i++) {
			double norm = Matrix2D.norm(eigenFacesData[i]);
			Log.d(LOG_TAG, "getNormalisedEigenFaces norm: " + norm);
			for (int j = 0; j < eigenFacesData[i].length; j++) {
				double v = eigenFacesData[i][j];
				eigenFacesData[i][j] = v / norm;
			}
		}
		return new Matrix2D(eigenFacesData);
	}

	private Matrix2D getNormalisedImagesData(Bitmap[] bufImgs) {
		int imageWidth = bufImgs[0].getWidth();
		int imageHeight = bufImgs[0].getHeight();

		int rows = bufImgs.length;
		int cols = imageWidth * imageHeight;

		int[][] data = new int[rows][cols];
		double[][] dataDouble = new double[rows][cols];

		for (int i = 0; i < rows; i++) {
			bufImgs[i].getPixels(data[i], 0, imageWidth, 0, 0, imageWidth,
					imageHeight);
		}

		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				dataDouble[i][j] = (double) (Color.blue(data[i][j]));
			}
		}

		Matrix2D imagesData = new Matrix2D(dataDouble);
		imagesData.normalise();
		return imagesData;
	}

	private void checkImageDimensions(List<String> filenames, Bitmap[] bufimgs)
			throws FaceRecError {
		int imgheight = 0;
		int imgwidth = 0;

		for (int i = 0; i < bufimgs.length; i++) {
			if (i == 0) {
				imgheight = bufimgs[i].getHeight();
				imgwidth = bufimgs[i].getWidth();
			}

			if ((imgheight != bufimgs[i].getHeight())
					|| (imgwidth != bufimgs[i].getWidth())) {
				String response = "all images should have same dimensions! "
						+ filenames.get(i) + " is of diff size";
				Log.e("FaceRec", response);
				throw new FaceRecError(response);
			}
		}
	}

	public Bitmap[] getGrayScaleImages(List<String> filenames)
			throws FaceRecError {
		Bitmap b = null;
		Bitmap[] bufimgs = new Bitmap[filenames.size()];

		Iterator<String> it = filenames.iterator();
		int i = 0;

		while (it.hasNext()) {
			String fn = it.next();
			File f = new File(fn);
			if (f.isFile()) {
				b = BitmapFactory.decodeFile(fn);

				if (b != null) {
					b = convertToGray(b);
					bufimgs[i++] = b;
				}
			}
		}
		return bufimgs;
	}

	private Bitmap convertToGray(Bitmap img) {
		int width = img.getWidth();
		int height = img.getHeight();

		Bitmap gray = Bitmap.createBitmap(width, height, Bitmap.Config.RGB_565);

		ColorMatrix cm = new ColorMatrix();
		cm.setSaturation(0);

		ColorMatrixColorFilter f = new ColorMatrixColorFilter(cm);

		Paint paint = new Paint();
		paint.setColorFilter(f);

		Canvas c = new Canvas(gray);
		c.drawBitmap(img, 0, 0, paint);

		return gray;
	}

	private void writeCache(String dir, FaceBundle cachedata) throws IOException {
		FileOutputStream fout = null;
		ObjectOutputStream fos = null;
		
		File file = new File(dir, "mycache.cache");
		if (file.exists())
			Log.d("FaceRec", "cache file exists, overwriting");
		
		fout = new FileOutputStream(file);
		fos = new ObjectOutputStream(fout);
		
		try {
			fos.writeObject(cachedata);
			Log.i("FaceRec", "wrote cache");
			
		} catch (Exception e) {
			e.printStackTrace();
			file.delete();
		}
		
		fout.close();
	}

	private FaceBundle getOldFacesBundle(String dir) throws IOException,
			ClassNotFoundException {
		File file = new File(dir, "mycache.cache");
		if (!file.exists()) {
			Log.e("FaceRec", "File " + file + " does not exist.");
		}
		
		FileInputStream fin = new FileInputStream(file);
		ObjectInputStream fo = new ObjectInputStream(fin);
		FaceBundle oldBundle = (FaceBundle) fo.readObject();
		fo.close();
		return oldBundle;
	}

	private void validateSelectedEigenFacesNumber(int selectedNumOfEigenFaces,
			List<String> newFileNames) throws FaceRecError {
		int numImgs = newFileNames.size();
		if (selectedNumOfEigenFaces <= 0 || selectedNumOfEigenFaces > numImgs) {
			Log.e("FaceRec", "incorrect number of selectedeigenfaces: "
					+ selectedNumOfEigenFaces + " used " + " allowed btw 0-" + numImgs);
			throw new FaceRecError(
					"incorrect number of selectedeigenfaces used..\n use a number between 0 and upto "
							+ (numImgs - 1));
		}
	}

	private List<String> getFileNames(String dir, String[] children) {
		java.util.List<String> imageFileNames = new java.util.ArrayList<String>();
		for (String i : children) {
			String fileName = dir + File.separator + i;
			imageFileNames.add(fileName);
		}
		Collections.sort(imageFileNames);
		return imageFileNames;
	}

	public List<String> parseDirectory(String directoryName, String extension)
			throws FaceRecError {
		final String ext = "." + extension;
		String[] children = null;
		File directory = new File(directoryName);

		if (directory.isDirectory()) {
			children = directory.list(new FilenameFilter() {
				public boolean accept(File f, String name) {
					return name.endsWith(ext);
				}
			});
		} else {
			throw new FaceRecError(directoryName + " is not a directory");
		}
		return getFileNames(directoryName, children);
	}

	public void checkCache(String dir, String extension,
			int selectedNumOfEigenFaces) throws FaceRecError, IOException {

		List<String> newFileNames = parseDirectory(dir, extension);
		FaceBundle oldBundle = null;

		try {
			validateSelectedEigenFacesNumber(selectedNumOfEigenFaces, newFileNames);
			oldBundle = getOldFacesBundle(dir);
			processCache(dir, newFileNames, oldBundle, selectedNumOfEigenFaces);
		} catch (FileNotFoundException e) {
			Log.i("FaceRec", "cache file not found");
			doCalculations(dir, newFileNames, selectedNumOfEigenFaces);
		} catch (Exception e) {
			e.printStackTrace();
			throw new FaceRecError(e.getMessage());
		}
	}

	private void processCache(String dir, List<String> newFileNames,
			FaceBundle oldBundle, int selectedNumOfEigenFaces) throws FaceRecError,
			IOException {
		List<String> oldFileNames = oldBundle.getImageFileNamesList();
		if (newFileNames.equals(oldFileNames)) {
			this.bundle = oldBundle;
			this.weights = calculateWeights(oldBundle, selectedNumOfEigenFaces);
		} else {
			Log.i("FaceRec", "folder contents changed");
			doCalculations(dir, newFileNames, selectedNumOfEigenFaces);
		}
	}

	private void createImageFromArray(String filename, double[] imgdata, int wd)
			throws IOException {
		Bitmap bmpGrayScale = Bitmap.createBitmap(wd, imgdata.length / wd,
				Bitmap.Config.RGB_565);

		double maxValue = Double.MIN_VALUE;
		double minValue = Double.MAX_VALUE;

		for (int i = 0; i < imgdata.length; i++) {
			maxValue = Math.max(maxValue, imgdata[i]);
			minValue = Math.min(minValue, imgdata[i]);
		}

		for (int j = 0; j < imgdata.length; j++) {
			imgdata[j] = ((imgdata[j] - minValue) * 255) / (maxValue - minValue);
		}

		int[] imgdataInt = new int[imgdata.length];

		// change double to int
		for (int index = 0; index < imgdata.length; index++)
			imgdataInt[index] = (int) imgdata[index];

		bmpGrayScale.setPixels(imgdataInt, 0, wd, 0, 0, wd, imgdata.length / wd);

		FileOutputStream out = new FileOutputStream(filename);
		bmpGrayScale.compress(Bitmap.CompressFormat.PNG, 100, out);
	}

	private String makeNewFolder(String fldr) {
		File folder = new File(fldr);

		if (folder.isDirectory()) {
			printError("folder:" + fldr + " exists");
			deleteContents(folder);

		} else {
			printError("no such folder as:" + fldr);
			boolean madeFolder = folder.mkdirs();
			if (!madeFolder) {
				printError("could not create folder :" + fldr);
			}

			// internal memory folder
			folder = this.androidContext.getDir(fldr, Context.MODE_PRIVATE);
			if (folder.exists()) {
				Log.d("FaceRec", "Folder created successfully in internal memory: "
						+ folder.getAbsolutePath());
			}
		}

		return folder.getAbsolutePath();
	}

	private void deleteContents(File f) {
		File[] files = f.listFiles();
		for (File i : files) {
			delete(i);
		}
	}

	private void delete(File i) {
		if (i.isFile()) {
			boolean deleted = i.delete();
			if (!deleted) {
				printError("file:" + i.getPath() + "could not be deleted");
			}
		}
	}

	public void reconstructFaces(int selectedeigenfaces) throws IOException {
		double[][] eigenfacesArray = bundle.getEigenFaces();
		Matrix2D egnfacesMatrix = new Matrix2D(eigenfacesArray);
		Matrix2D egnfacesSubMatrix = egnfacesMatrix
				.getSubMatrix(selectedeigenfaces);
		double[] eigenvalues = bundle.getEigenValues();
		Matrix2D eigenvalsMatrix = new Matrix2D(eigenvalues, 1);
		Matrix2D eigenvalsSubMatrix = eigenvalsMatrix.transpose().getSubMatrix(
				selectedeigenfaces);

		/*
		 * the term 'phi' is used to denote mean subtracted reconstructed imagedata
		 * since that term appears in T&P doc.The term 'xnew' denotes
		 * phi+average_image
		 */
		double[][] phi = getPhiData(egnfacesSubMatrix, eigenvalsSubMatrix);
		double[][] xnew = addAverageFaceData(phi);
		String reconFolderName = getAbsoluteFolderPath("reconfaces");
		String ext = ".png";
		reconstructPhiImages(phi, reconFolderName, ext);
		reconstructOriginalImages(xnew, reconFolderName, ext);
		Log.i("FaceRec", "reconstruction over");
	}

	private double[][] getPhiData(Matrix2D egnfacesSubMatrix,
			Matrix2D eigenvalsSubMatrix) {
		double[] evalsub = eigenvalsSubMatrix.flatten();
		Matrix2D tempEvalsMat = new Matrix2D(weights.length, evalsub.length);
		tempEvalsMat.replaceRowsWithArray(evalsub);
		Matrix2D tempmat = new Matrix2D(weights);
		tempmat.multiplyElementWise(tempEvalsMat);
		Matrix2D phinewmat = tempmat.multiply(egnfacesSubMatrix);
		double[][] phi = phinewmat.toArray();
		return phi;
	}

	private double[][] addAverageFaceData(double[][] phi) {
		double[][] xnew = new double[phi.length][phi[0].length];
		double[] avgface = bundle.getAvgFace();
		for (int i = 0; i < phi.length; i++) {
			for (int j = 0; j < phi[i].length; j++) {
				xnew[i][j] = phi[i][j] + avgface[j];
			}
		}
		return xnew;
	}

	private void reconstructOriginalImages(double[][] xnew,
			String reconFolderName, String ext) throws IOException {
		String prefix;
		prefix = "xnew";
		int imgwidth = bundle.getImageWidth();
		for (int i = 0; i < xnew.length; i++) {
			double[] xnewdata = xnew[i];
			String filename = reconFolderName + File.separator + prefix + i + ext;
			createImageFromArray(filename, xnewdata, imgwidth);
		}
	}

	private void reconstructPhiImages(double[][] phi, String reconFolderName,
			String ext) throws IOException {
		int imgwidth = bundle.getImageWidth();
		reconFolderName = makeNewFolder(reconFolderName);
		String prefix = "phi";
		for (int i = 0; i < phi.length; i++) {
			double[] phidata = phi[i];
			String filename = reconFolderName + File.separator + prefix + i + ext;
			createImageFromArray(filename, phidata, imgwidth);
		}
	}

	private int getNumofFacesVal(String numofFaces) throws NumberFormatException {
		return Integer.valueOf(numofFaces);
	}

	private double getThresholdVal(String threshold) throws NumberFormatException {
		return Double.valueOf(threshold);
	}

	private void validateSelectedImageFileName(String faceimagename)
			throws FaceRecError {
		if ((faceimagename.length() == 0) || (!isImage(faceimagename))) {
			throw new FaceRecError("select an imagefile");
		}
	}

	private boolean isImage(String imagefilename) {
		return Utils.isImageFile(imagefilename);
	}

	private void validateSelectedFolderName(String foldername)
			throws FaceRecError {
		if (foldername.length() == 0) {
			throw new FaceRecError("select a folder");
		}
	}

	public MatchResult processSelections(String faceImageName, String directory,
			String numofFaces, String threshold) {

		MatchResult result = null;
		int numFaces = 0;
		double thresholdVal = 0.0;

		try {
			Log.d(LOG_TAG, "processSelections+");
			
			Log.d(LOG_TAG, "processSelections: " + "validation");
			validateSelectedImageFileName(faceImageName);
			validateSelectedFolderName(directory);

			
			numFaces = getNumofFacesVal(numofFaces);
			thresholdVal = getThresholdVal(threshold);
			String extension = getFileExtension(faceImageName);
			Log.d(LOG_TAG, "processSelections: " + "numFaces: " + numFaces);
			Log.d(LOG_TAG, "processSelections: " + "thresholdVal: " + thresholdVal);
			Log.d(LOG_TAG, "processSelections: " + "extension: " + extension);

			
			checkCache(directory, extension, numFaces);
			reconstructFaces(numFaces);

			result = findMatchResult(faceImageName, numFaces, thresholdVal);
			Log.d(LOG_TAG, "processSelections: " + "result: " + result.getMatchSuccess());
			Log.d(LOG_TAG, "processSelections: " + "result: " + result.getMatchDistance());
			Log.d(LOG_TAG, "processSelections: " + "result: " + result.getMatchMessage());
			Log.d(LOG_TAG, "processSelections: " + "result: " + result.getMatchFileName());

		} catch (Exception e) {
			Log.e("FaceRec", "Exception thrown..");
			e.printStackTrace();
			result = new MatchResult(false, null, Double.NaN, e.getMessage());
		}
		return result;
	}

	public static String getFileExtension(String filename) {
		String ext = "";
		int i = filename.lastIndexOf('.');
		if (i > 0 && i < filename.length() - 1) {
			ext = filename.substring(i + 1).toLowerCase();
		}
		return ext;
	}

	public static void printError(String msg) {
		Log.e("FaceRec", msg);
	}

	public static void debug(String msg) {
		Log.d("FaceRec", msg);
	}
	
	private String getAbsoluteFolderPath(String folder) {
		return this.androidContext.getDir(folder, Context.MODE_PRIVATE).getAbsolutePath();
	}
}
