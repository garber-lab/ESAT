package broad.pda.enrichment;

import java.io.File;
import java.util.TreeMap;
import java.util.SortedMap;

import org.apache.commons.io.FilenameUtils;

public class EnrichmentUtils {
	
	/**
	 * Recursively visits directories to find genomic annotation files.
	 * @param dir		Base directory
	 * @param fileExt	List of valid genomic annotation file extensions
	 * @return			Map of genomic annotation names and the corresponding files.
	 */
	public static SortedMap<String, File> findAnnotationFiles(File dir, String[] fileExt) {
		SortedMap<String,File> map = new TreeMap<String,File>();
		findAnnotationFiles(dir, fileExt, "", map);
		return map;
	}
	
	public static SortedMap<String, File> findAnnotationFiles(File dir, String[] fileExt, String basename) {
		SortedMap<String,File> map = new TreeMap<String,File>();
		findAnnotationFiles(dir, fileExt, basename, map);
		return map;
	}
	
	public static SortedMap<String, File> findAnnotationFiles(File dir, String[] fileExt, String basename, SortedMap<String,File> map) {
		if (dir.isFile()) {
			// Then we've reached end of the recursion
			for (int i = 0; i < fileExt.length; i++) {
				if (FilenameUtils.getExtension(dir.getName()).equals(fileExt[i])) {
					String name = basename + FilenameUtils.removeExtension(dir.getName());
					map.put(name, dir);
					break;
				}
			}
		} else if (dir.isDirectory()) {
			File[] files = dir.listFiles();
			for (File file : files) {
				findAnnotationFiles(file, fileExt, basename + FilenameUtils.getName(dir.getAbsolutePath()) + ".", map);
			}
		}
		
		return map;
	}
	
	
}
