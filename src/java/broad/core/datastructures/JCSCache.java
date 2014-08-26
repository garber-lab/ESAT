package broad.core.datastructures;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Random;
import java.util.TreeMap;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.feature.GenomeWindow;
import nextgen.core.feature.Window;
import nextgen.core.model.JCSAlignmentModel.AlignmentCount;
import nextgen.core.model.JCSAlignmentModel.FilteredIterator;
import nextgen.core.model.JCSAlignmentModel;
import nextgen.core.readers.PairedEndReader;

import org.apache.jcs.JCS;
import org.apache.jcs.access.exception.CacheException;
import org.apache.jcs.engine.control.CompositeCacheManager;
import org.apache.log4j.Logger;
import org.apache.log4j.Level;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;

public class JCSCache {
	
	static Logger logger = Logger.getLogger(JCSCache.class.getName());
	/** Default directory where cache files will be spooled: ~/jcs **/
	private final String defaultDir = System.getProperty("java.io.tmpdir") + "/.jcs";
	
//	private static final String DEFAULT_GROUP = "Default_JCS_Group";
	private static final String DEFAULT_NAME = "Default_JCS_Name";
	private static final int KEY_LENGTH = 9;

	static Collection<Integer> keys;
	private final Long maxLifeSeconds=7200L;
	/** Maximum number of in-memory instances before sending items to disk. Default is 50,000. */
	private Long defaultCapacity = 200000L;
	IntervalTree<Integer> keyTree;
	PairedEndReader reader;
	String cacheChr = null;
	int cacheSize;
	int cacheStart=0;
	int cacheEnd = 0;
	boolean fullyContained = false;
	JCSAlignmentModel model;
	//String groupName;
	String cacheName;
//	private static JCS treeCache;
	
	/** Holds a list of already defined caches to help ensure uniqueness. */
	private static List<String> cacheNames = new ArrayList<String>();
	
	/** A thread that will ensure that all of these caches will be disposed of during shutdown. */ 
	private static Thread shutdownThread = new Thread() {
		/** Run the shutdown hook for disposing of all caches. */
		@Override 
		public void run() {
			for (String cacheName : cacheNames) {
				try {
		              logger.info("Shutting down " + cacheName + " cache.");
		              JCS.getInstance(cacheName).clear();
		              JCS.getInstance(cacheName).dispose();
		            }
	            catch (Exception e) {
	              String msg = "Failure to clear cache " + cacheName + ":" + e.getMessage();
	              logger.info(msg);
	            }
	          }
	        }
	      };
	/**
	 * Initializes the cache
	 * @param reader
	 */
	public JCSCache(PairedEndReader reader,int cacheSize,JCSAlignmentModel m){
		synchronized (JCSCache.class) {
			this.reader=reader;
			keyTree = new IntervalTree<Integer>();
			
			//SET THE DEFAULT CAPACITY
			//EACH RECORD TAKES APPROX 348728 BYTES IN MEMORY
			//defaultCapacity = (long)(Runtime.getRuntime().maxMemory()/(2.0*175000));
			Logger.getLogger("org.apache.jcs").setLevel(Level.OFF);
			cacheName=JCSCache.DEFAULT_NAME+"_"+System.currentTimeMillis();
			CompositeCacheManager ccm = CompositeCacheManager.getUnconfiguredInstance();
			ccm.configure(initJcsProps(cacheName));
			cacheNames.add(cacheName);
//			Runtime.getRuntime().addShutdownHook(JCSCache.shutdownThread);
			
			this.cacheSize = cacheSize;
			this.model=m;
			mkdirs(defaultDir);

			keys = new HashSet<Integer>();
			logger.info("Indexed disk cache written to: "+defaultDir);
			logger.info("Holding "+defaultCapacity+" in In-memory cache");
/*
			try{
				treeCache = JCS.getInstance(cacheName);
			}catch (Exception e)
	        {
	            // Handle cache region initialization failure
	        }*/
			
		}
	}
	
	/**
	 * Generates a Long random number of length "length". We will use 10 digits
	 * @param length
	 * @return
	 */
	public static Integer generateRandom(int length) {
	    Random random = new Random();
	    char[] digits = new char[length];
	    digits[0] = (char) (random.nextInt(9) + '1');
	    for (int i = 1; i < length; i++) {
	        digits[i] = (char) (random.nextInt(10) + '0');
	    }
	    if(keys.contains(Integer.parseInt(new String(digits)))){
	    	return generateRandom(length);
	    }
	    return Integer.parseInt(new String(digits));
	}
	/**
	 * 
	 * @return
	 */
	private Properties initJcsProps(String cacheName) {
		
		String reg = "jcs.region." + cacheName;
	    String regCacheAtt = reg + ".cacheattributes";
	    String regEleAtt = reg + ".elementattributes";
	    //String aux = "jcs.auxiliary.DC-" + cacheName;
	    String aux = "jcs.auxiliary.DC";
	    String auxAtt = aux + ".attributes";
	    String elementAtt = "jcs.default.elementattributes";
	    String memName = "org.apache.jcs.engine.memory.lru.LRUMemoryCache";
	    String diskAttName = "org.apache.jcs.auxiliary.disk.indexed.IndexedDiskCacheAttributes";
	    Properties props = new Properties();
	    props.setProperty("jcs.default", "DC");
	    //props.setProperty("jcs.default", "");
	    //props.setProperty("jcs.default.cacheattributes.MaxObjects", defaultCapacity.toString());
	    //props.setProperty("jcs.default.elementattributes.MaxLifeSeconds", "15000");
	    props.setProperty("jcs.default.cacheattributes.DiskUsagePatternName", "SWAP");

	    props.setProperty(elementAtt + ".IsEternal", "false");
	    props.setProperty(elementAtt + ".IsSpool", "true");
	    props.setProperty(elementAtt + ".IsRemote", "false");
	    props.setProperty(elementAtt + ".IsLateral", "false");

	    //props.setProperty(reg, "DC-" + cacheName);
	    props.setProperty(reg, "DC");
	    //props.setProperty(reg, "" + cacheName);
	    props.setProperty(regCacheAtt, "org.apache.jcs.engine.CompositeCacheAttributes");
	    props.setProperty(regCacheAtt + ".MaxObjects", defaultCapacity.toString());
	    props.setProperty(regCacheAtt + ".MemoryCacheName", memName);
	    props.setProperty(regCacheAtt + ".UseMemoryShrinker", "true");
	    props.setProperty(regCacheAtt + ".MaxMemoryIdleTimeSeconds", "15000");
	    props.setProperty(regCacheAtt + ".ShrinkerIntervalSeconds", "7200");
	    props.setProperty(regCacheAtt + ".DiskUsagePatternName", "SWAP");
	   // props.setProperty(regCacheAtt + ".MaxSpoolPerRun", "5000000");
	    props.setProperty(regEleAtt, "org.apache.jcs.engine.ElementAttributes");
	    props.setProperty(regEleAtt + ".IsEternal", "false");
	    props.setProperty(regEleAtt + ".IsSpool", "true");
	    props.setProperty(regEleAtt + ".IsRemote", "false");
	    props.setProperty(regEleAtt + ".IsLateral", "false");
	    props.setProperty(regEleAtt + ".MaxLifeSeconds", maxLifeSeconds.toString());
	    props.setProperty(aux, "org.apache.jcs.auxiliary.disk.indexed.IndexedDiskCacheFactory");
	    props.setProperty(auxAtt, diskAttName);
	    props.setProperty(auxAtt + ".DiskPath", defaultDir);
	    props.setProperty(auxAtt + ".maxKeySize", "100000000");
	    props.setProperty(auxAtt + ".OptimizeAtRemoveCount", "100000");
	    props.setProperty(auxAtt + ".ClearDiskOnStartup", "true");
	    props.setProperty(auxAtt + ".OptimizeOnShutdown", "true");
	    props.setProperty(auxAtt + ".MaxPurgatorySize", "10000");
	    //writeProperties(props);
	    return props;
	}
	
	  /**
	   * Returns a string containing all current properties of the JCSCache in alphabetical order.
	   * For logging and debugging. 
	   * @param properties The properties of interest.
	   */
	  private void writeProperties(Properties properties) {
	    String cr = System.getProperty("line.separator"); 
	    String eq = "=";
	    // Adding them to a treemap has the effect of alphabetizing them. 
	    TreeMap<String, String> alphaProps = new TreeMap<String, String>();
	    for (Map.Entry<Object, Object> entry : properties.entrySet()) {
	      String propName = (String)entry.getKey();
	      String propValue = (String)entry.getValue();
	      alphaProps.put(propName, propValue);
	    }
	    StringBuffer buff = new StringBuffer(25);
	    for (String key : alphaProps.keySet()) {
	      buff.append(key).append(eq).append(properties.get(key)).append(cr);
	    }
	    // Now print the properties to a file. 
	    BufferedWriter f;
	    try {
	      long timestamp = new Date().getTime();
	      f = new BufferedWriter(new FileWriter("jcs." + timestamp + ".properties"));
	      f.write(buff.toString());
	      f.close();
	    }
	    catch (Exception e) {
	      e.printStackTrace();
	    }
	  }


	/**
	 * Instantiates the directory indicated by dir, and throws a Runtime exception if unsuccessful.
	 * 
	 * @param dir The directory where the cache files will go.
	 */  
	private void mkdirs(String dir) {
		File path = new File(dir);
		boolean done = path.mkdirs();
		if (!done && !path.exists()) {
			throw new RuntimeException("mkdirs() failed");
		}
	}
	
	/**
	 * Return an iteratory over an annotation window
	 * @param window
	 * @param fullyContained
	 * @return
	 * @throws CacheException 
	 */
	private CloseableIterator<AlignmentCount> query(Annotation window, boolean fullyContained) throws CacheException {
		//if larger than the cache size then just return the query directly
		if(window.getSize()>this.cacheSize){
			//logger.info("Get reads for the entire window of size "+window.getSize()+" for "+window.toUCSC());
			return getReads(window, fullyContained);
		}
		//else if doesnt contain the window then update cache and query again
		else if (!contains(window) || this.fullyContained != fullyContained) {
			//logger.info("this.fullyContained = " + this.fullyContained + "  fullyContained = " + fullyContained);
			//logger.info("Updating cache for "+window.getSize()+" for "+window.toUCSC());
			updateCache(window.getReferenceName(), window.getStart(), window.getEnd(), fullyContained);
		}
		//pull reads from cache
		return getReadsFromCache(window);
	}
	
	public FilteredIterator query(Annotation window, boolean fullyContained, CoordinateSpace cs)  {
		try {
			CloseableIterator<AlignmentCount> iter = query(window, fullyContained);
//			logger.info("In query for window:"+window.toBED());
			FilteredIterator filtered=model.new FilteredIterator(iter, window, cs, fullyContained);
			return filtered;
		} catch (CacheException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			
//			FilteredIterator filtered = model.new FilteredIterator(new CloseableIterator(Collections.<AlignmentCount>emptyList().iterator()));
			//TODO: FIX
			return null;
		}
		
	}
	
	private CloseableIterator<AlignmentCount> getReadsFromCache(Annotation window) {
		// TODO:  Can't we just call cachedTree.overlappersValueIterator(window.getStart(), window.getEnd()) and get rid of the NodeIterator?
		//  it looks like overlappingValueIterator can handle multiple values per node
		
		//TODO
		//CODE return over alignment count
		
		return new JCSNodeIterator(this.keyTree.overlappers(window.getStart(), window.getEnd()));
	}
	
	public class JCSNodeIterator implements CloseableIterator<AlignmentCount>{
		Iterator<Node<Integer>> iter;

		JCSNodeIterator (Iterator<Node<Integer>> overlappers) {
			iter=overlappers;
		}

		@Override
		public boolean hasNext() {
			return iter.hasNext();
		}

		@Override
		public AlignmentCount next() {
			Node<Integer> keyNode=iter.next();
			Collection<Alignment> reads = new HashSet<Alignment>();
			boolean flag=true;
			Alignment value=null;
			for(Integer key:keyNode.getContainedValues()){
			try{
//				for(Object key1:JCS.getInstance(cacheName).getGroupKeys(groupName)){
					//logger.info((String)key1);
//				}
//				System.out.println("Attempt to GET "+key);
				//Alignment read = (Alignment)treeCache.getFromGroup((Object)key,groupName);
				Alignment read = (Alignment)JCS.getInstance(cacheName).get((Object)key);

				if(read!=null){
					if(read.getHeader()==null)
						read.setHeader(model.getHeader());
					reads.add(read);
					//logger.info("Read is read");
				}
				else{
				
/*					if(keys.contains(key)){
						if(JCS.getInstance(cacheName).get(key)==null){
							logger.info("Read is NULL BUT CONTAINS the key");
						}
						else{
							if(keys.contains(key)){
								logger.info("Read is NULL BUT DOES NOT contain the key but KEYS contains the KEY Cache "+cacheStart+":"+cacheEnd+" VS "+key);
							}
							else
								logger.info("Read is NULL and its OKAY.");
						}
					}*/
				}				
				if(flag && read!=null){
					flag=false;
					value = read;
				}
			}catch(CacheException e){
				try {
					JCS.getInstance(cacheName).clear();
				} catch (CacheException e1) {
					logger.error(e.getMessage());
				}
				logger.error(e.getMessage());
			}
			}
			return model.new AlignmentCount(value,reads);
		}

		@Override
		public void remove() {}

		@Override
		public void close() {}
	}

	/**
	 * Update the cache to have these new positions
	 * @param chr
	 * @param start
	 * @param end
	 * @throws CacheException 
	 */
	private void updateCache(String chr, int start, int end, boolean fullyContained) throws CacheException {
		int newStart=start;
		int newEnd=end;

		logger.debug("Updating cache: " + chr + ":" + start + "-" + end + ". CacheSize=" + cacheSize);
		// if window is larger than cache size 
		//@skadri TODO: Isn't this checked in query() already?
		// (this will happen in TranscriptomeSpace if a transcript is longer than the cache size)
		if ((end-start) > this.cacheSize) {
			newStart=start;
			newEnd=end;
		} else {
			if (this.cacheChr == null || chr.equalsIgnoreCase(this.cacheChr)) { 
				if (end > this.cacheEnd) {
					// Set the cache to start and the beginning of the requested sequence
					newStart = start;
					newEnd = start + this.cacheSize;
				} else if (start < this.cacheStart) {
					// Maybe we're scanning backwards?  So we'll fix the cache to the end of the window
					newEnd = end;
					newStart = end - this.cacheSize;

//					throw new IllegalStateException("Cache DEBUG");
				}
			}
		}

		this.cacheChr=chr;
		this.cacheStart=newStart;
		this.cacheEnd=newEnd;
		this.fullyContained = fullyContained;

		Window update=new GenomeWindow(this.cacheChr, this.cacheStart, this.cacheEnd);

		this.keyTree=getIntervalTree(update, fullyContained);
		
		logger.debug("Cache updated to: " + cacheChr + ":" + cacheStart + "-" + cacheEnd);
	}
	
	/**
	 * Returns an interval tree of reads over the specified window
	 * @param w
	 * @param fullyContained
	 * @return
	 * @throws CacheException 
	 */
	private IntervalTree<Integer> getIntervalTree(Window w, boolean fullyContained) throws CacheException {
		
		//logger.info("Get interval tree");
		//Assume update cache will not fail
		//Set at 2 million reads
		//CLEAN THE CACHE
/*		if(groupName!=null){
			for(Object key:treeCache.getGroupKeys(groupName)){
				treeCache.remove(key, groupName);
			}
		}	*/
/*		if(cacheName!=null)
			dispose(cacheName);*/
		
		JCS.getInstance(cacheName).clear();
		JCS.getInstance(cacheName).dispose();
		
/*		for(Integer key:keys){
			Alignment a = (Alignment)JCS.getInstance(cacheName).get((Object)key);
			if(a!=null){
				logger.info("Cache still contains object after clearing.");
			}
		}*/
/*		if(!keys.isEmpty()){
			for(String k:keys){
				try {
					treeCache.remove((Object)k);
				} catch (CacheException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}*/
		
	 	IntervalTree<Integer> tree=new IntervalTree<Integer>();
	 	keys = new HashSet<Integer>();
//	 	groupName = this.DEFAULT_GROUP+w.getChr()+w.getStart();
//	 	cacheName = DEFAULT_NAME+System.currentTimeMillis();
//	 	cacheNames.add(cacheName);
	 	//Initialize the cache
//		CompositeCacheManager ccm = CompositeCacheManager.getUnconfiguredInstance();
//		ccm.configure(initJcsProps(cacheName));
	 	
		CloseableIterator<AlignmentCount> iterReadsOverlappingRegion=getReads(w, fullyContained);
		while(iterReadsOverlappingRegion.hasNext()){
/*			final int unit = 1024; // kB
			long prevusedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		    long availableMemory = Runtime.getRuntime().maxMemory() - prevusedMemory;
		    logger.info("Memory: free="+(Runtime.getRuntime().freeMemory()/unit)+" total="+(Runtime.getRuntime().totalMemory()/unit)+" max="+(Runtime.getRuntime().maxMemory()/unit+" used="+prevusedMemory/unit+" available="+availableMemory/unit));
*/
			Alignment record=iterReadsOverlappingRegion.next().getRead();
/*		    long usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
		    availableMemory = Runtime.getRuntime().maxMemory() - usedMemory;
		    logger.info("Memory: free="+(Runtime.getRuntime().freeMemory()/unit)+" total="+(Runtime.getRuntime().totalMemory()/unit)+" max="+(Runtime.getRuntime().maxMemory()/unit+" used="+usedMemory/unit+" available="+availableMemory/unit));
		    logger.info("Object adding used :"+(usedMemory-prevusedMemory)/unit);
*/
			if (model.isValid(record)) {
				//logger.info("Attempting to ADD "+record.getName());
				//String key = record.getReadName()+"-"+record.getAlignmentStart()+":"+record.getAlignmentEnd();
				//logger.info(record.getName());
				
				Integer key = generateRandom(KEY_LENGTH);

				//				long start = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
				JCS.getInstance(cacheName).put(key,record);
			    keys.add(key);    
//			    long end = Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
//				if((end-start)!=0){
//					logger.info(end-start);
//				}

				tree.put(record.getAlignmentStart(), record.getAlignmentEnd(), key);
			}	
			/*Node<Alignment> node=tree.find(record.getAlignmentStart(), record.getAlignmentEnd());			
			if(node!=null){node.incrementCount();}
			else{tree.put(record.getAlignmentStart(), record.getAlignmentEnd(), record);}*/
		}			
		iterReadsOverlappingRegion.close();
		/*		logger.info("Checking the tree");
		for(String k:keys){
			
			
			try {
				Alignment read = (Alignment)JCS.getInstance(cacheName).get((Object)k);
				if(read==null){
			    	  logger.info("Read just added is returning null "+k);
			      }
			      else{
			    	  logger.info("Read is OK "+k);
			      }
			} catch (CacheException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		      
		}*/
		return tree;
	}
	
	
	public CloseableIterator<AlignmentCount> getReads(Annotation w, boolean fullyContained){
		return model.new WrapAlignmentCountIterator(this.reader.query(w, fullyContained));
	}
	
	
	public CloseableIterator<AlignmentCount> getReads(){
		return model.new FilteredIterator(model.new WrapAlignmentCountIterator(reader.iterator()));
	}
	
	/**
	 * Shuts down the specified cache, and removes it from the list of active caches so it can be
	 * created again.
	 * 
	 * @param cacheName The name of the cache to dispose of.
	 */
	public static void dispose(String cacheName) {
        try {
          cacheNames.remove(cacheName);
          JCS.getInstance(cacheName).clear();
        }
        catch (CacheException e) {
          String msg = "Failure to clear cache " + cacheName + ":" + e.getMessage();
          logger.error(msg);
        }
      }     
	
	/**
	 * Does the current cache have the window?
	 * @param region The window to determine whether its in the cache
	 */
	private boolean contains(Annotation region) {
		if (this.cacheChr==null) {
			//Not yet initialized so cant be contained in cache
			return false;
		}
		
		if (this.cacheChr.equalsIgnoreCase(region.getReferenceName())) {
			if (this.cacheStart <= region.getStart() && this.cacheEnd >= region.getEnd()) { return true; }
		}
		
		return false;
	}
	
}
