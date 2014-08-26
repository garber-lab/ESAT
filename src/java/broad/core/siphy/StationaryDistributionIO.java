package broad.core.siphy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.math.EmpiricalDistribution;
import broad.core.siphy.EvolutionaryModel.PiFit;

public class StationaryDistributionIO {
	IntervalTree<PiFit> fits;
	//ArrayList<PiFit> fits;
	
	public StationaryDistributionIO () {
		super();
		//fits = new ArrayList<PiFit>();
		fits = new IntervalTree<PiFit>();
	}
	
	public void loadIntoDistribution(File source, EmpiricalDistribution ed, boolean shortFormat) throws IOException {
		FileInputStream fis = new FileInputStream(source);
		try {
			load(fis);
		} finally {
			fis.close();
		}
	}
	public void load(InputStream is) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
			}
			String [] lineInfo = line.trim().split("\t");
			PiFit fit = new PiFit(lineInfo);
			fits.put(fit.getPosition(),fit.getPosition() + 1,fit);
		}
	}
	
	public void load(String filePath) throws IOException {
		FileInputStream fis = new FileInputStream(filePath);
		load(fis);
		fis.close();
	}
	

	public void shift(int amountToShift) {
		Iterator<PiFit> it = fits.valueIterator();
		IntervalTree<PiFit> newTree = new IntervalTree<PiFit>(); //It was too slow (and actually not working) to remove and add to the old tree
 		while(it.hasNext()) {
			PiFit fit = it.next();
			fit.setPosition(fit.getPosition() + amountToShift);
			/*
			fits.remove(fit.getPosition(), fit.getPosition() + 1);
			fit.setPosition(fit.getPosition() + amountToShift);
			fits.put(fit.getPosition(), fit.getPosition() + 1, fit);
			*/
			newTree.put(fit.getPosition(), fit.getPosition() + 1, fit);
		}
 		fits = newTree;
	}
	
	public List<PiFit> getFits() {
		List<PiFit> fitList = new ArrayList<PiFit>(fits.size());
		Iterator<PiFit> fitIt = fits.valueIterator();
		while(fitIt.hasNext()) {
			fitList.add(fitIt.next());
		}
		return fitList;
	}
	
	public Iterator<PiFit> getOverlappers(LightweightGenomicAnnotation annotation) {
		return new IntervalTree.ValuesIterator<PiFit>(fits.overlappers(annotation.getStart(), annotation.getEnd()));
	}

}
