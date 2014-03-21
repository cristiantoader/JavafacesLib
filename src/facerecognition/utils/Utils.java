package facerecognition.utils;

import android.graphics.BitmapFactory;

public class Utils {
	public static boolean isImageFile(String filename) {
		BitmapFactory.Options options = null;
		
		if (filename == null || filename.length() == 0) {
			return false;
		}
		
		options = new BitmapFactory.Options();
		options.inJustDecodeBounds = true;
		BitmapFactory.decodeFile(filename);
		
		return options.outWidth != -1 && options.outHeight != -1;
	}
}
