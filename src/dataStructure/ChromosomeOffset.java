package dataStructure;

import java.util.HashMap;

/**
 * 
 * 
 */
public enum ChromosomeOffset {

	/**
	 * 
	 */
	CHR_1((short) 1, "chr1", 0 + 100000),
	CHR_2((short) 2, "chr2", 254235634 + CHR_1.getOffset() + 100000),
	CHR_3((short) 3, "chr3", 248063361 + CHR_2.getOffset() + 100000),
	CHR_4((short) 4, "chr4", 201982879 + CHR_3.getOffset() + 100000),
	CHR_5((short) 5, "chr5", 194977362 + CHR_4.getOffset() + 100000),
	CHR_6((short) 6, "chr6", 184533566 + CHR_5.getOffset() + 100000),
	CHR_7((short) 7, "chr7", 174537369 + CHR_6.getOffset() + 100000),
	CHR_8((short) 8, "chr8", 162321437 + CHR_7.getOffset() + 100000),
	CHR_9((short) 9, "chr9", 149291303 + CHR_8.getOffset() + 100000),
	CHR_10((short) 10, "chr10", 144037700 + CHR_9.getOffset() + 100000),
	CHR_11((short) 11, "chr11", 138245442 + CHR_10.getOffset() + 100000),
	CHR_12((short) 12, "chr12", 137706647 + CHR_11.getOffset() + 100000),
	CHR_13((short) 13, "chr13", 136528933 + CHR_12.getOffset() + 100000),
	CHR_14((short) 14, "chr14", 117473276 + CHR_13.getOffset() + 100000),
	CHR_15((short) 15, "chr15", 109496531 + CHR_14.getOffset() + 100000),
	CHR_16((short) 16, "chr16", 104582020 + CHR_15.getOffset() + 100000),
	CHR_17((short) 17, "chr17", 92161849 + CHR_16.getOffset() + 100000),
	CHR_18((short) 18, "chr18", 82819115 + CHR_17.getOffset() + 100000),
	CHR_19((short) 19, "chr19", 79638793 + CHR_18.getOffset() + 100000),
	CHR_20((short) 20, "chr20", 60311563 + CHR_19.getOffset() + 100000),
	CHR_21((short) 21, "chr21", 64286031 + CHR_20.getOffset() + 100000),
	CHR_22((short) 22, "chr22", 49092493 + CHR_21.getOffset() + 100000),
	CHR_23((short) 23, "chrX", 52330658 + CHR_22.getOffset() + 100000),
	CHR_24((short) 24, "chrY", 158375972 + CHR_23.getOffset() + 100000),
	CHR_25((short) 25, "chrM", 60561038 + CHR_24.getOffset() + 100000);
	

	/**
	 * Contains the chromosome number again
	 */
	@SuppressWarnings("unused")
	private final short num;
	/**
	 * contains the offset / start point of the chromosome relative to a
	 * concatenated genome
	 */
	private final long offset;
	
	private final String name;
	
	private static HashMap<String, ChromosomeOffset> name_to_offset;
	private static HashMap<Short, ChromosomeOffset> num_to_offset; 
	
	static {
		name_to_offset = new HashMap<String, ChromosomeOffset>();
		num_to_offset = new HashMap<Short, ChromosomeOffset>();
		
		for (ChromosomeOffset co : ChromosomeOffset.values()){
			name_to_offset.put(co.name.toLowerCase(), co);
			num_to_offset.put(co.num, co);
		}
		
	}

	/**
	 * @param num
	 *            chromosome number
	 * @param offset
	 *            relative genome position
	 */
	ChromosomeOffset(short num, String name, long offset) {
		this.num = num;
		this.offset = offset;
		this.name = name;
	}

	
	/**
	 * @return the relative genome position of the chromosome
	 */
	public long getOffset() {
		return this.offset;
	}
	
	/**
	 * @return the number of the chromosome according to the cmap file
	 */
	public short getChromosomeNumber() {
		return this.num;
	}
	
	/**
	 * @return the original name of the chromosome
	 */
	public String getChromosomeName() {
		return this.name;
	}
	
	/**
	 * Search the Offsets by chr number
	 * @param num short number of chr
	 * @return ChromosomeOffset
	 */
	public static ChromosomeOffset getChromosomeOffsetbyNumber(short num){
		return num_to_offset.get(num);
	}

	/**
	 * Method to get the chromosome number of a given chromosome name
	 * @param name String that represents the chromosome name
	 * @return the short number of the chromosome
	 */
	public static short chromosomeNumber(String name) {
		name = name.toLowerCase();
		if(name_to_offset.containsKey(name)){
			return name_to_offset.get(name).num;
		} else if(name.startsWith("chr")){ 
			return num_to_offset.get(Short.parseShort(name.substring(3))).num;
		} else if (! name.startsWith("chr")) {
			return name_to_offset.get("chr"+name).num;		
		} else {
			return num_to_offset.get(Short.parseShort(name)).num;
		}
	}
	
	/**
	 * @param i
	 *            chromosome number
	 * @return relative genome position of i-th chromosome
	 */
	public static long offset(short i) {
		return num_to_offset.get(i).offset;
	
	}
}