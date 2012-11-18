package dataStructure;

public class GATK_SNV_entry {
	private String chr;
	private Long genomPosition;
	private Long absGenomPosition;
	private Short chr_nr;
	private String refAllel;
	private String altAllel;
	private String dbSNP_entry;
	private Double qual_SNV;
	private String filterInfo;
	private String gatkInfo;
	private String gatkFormat;
	private String gatkGenotype;

	public GATK_SNV_entry(String chr, Long genomPosition,
			Long absGenomPosition, Short chr_nr) {
		this.chr = chr;
		this.genomPosition = genomPosition;
		this.absGenomPosition = absGenomPosition;
		this.chr_nr = chr_nr;
	}
	
	public GATK_SNV_entry(String chr, Long genomPosition,
			Long absGenomPosition, Short chr_nr, String dbSNP_entry, String refAllel,
			String altAllel, Double qual_SNV,
			String filterInfo, String gatkInfo, String gatkFormat,
			String gatkGenotype) {
		this(chr, genomPosition, absGenomPosition,chr_nr);
		this.dbSNP_entry = dbSNP_entry;
		this.refAllel = refAllel;
		this.altAllel = altAllel;
		this.qual_SNV = qual_SNV;
		this.filterInfo = filterInfo;
		this.gatkInfo = gatkInfo;
		this.gatkFormat = gatkFormat;
		this.gatkGenotype = gatkGenotype;
	}


	public String getChr() {
		return chr;
	}

	public Long getGenomPosition() {
		return genomPosition;
	}

	public Long getAbsGenomPosition() {
		return absGenomPosition;
	}

	public Short getChr_nr() {
		return chr_nr;
	}

	public String getRefAllel() {
		return refAllel;
	}

	public String getAltAllel() {
		return altAllel;
	}

	public String getDbSNP_entry() {
		return dbSNP_entry;
	}

	public Double getQual_SNV() {
		return qual_SNV;
	}

	public String getFilterInfo() {
		return filterInfo;
	}

	public String getGatkInfo() {
		return gatkInfo;
	}

	public String getGatkFormat() {
		return gatkFormat;
	}

	public String getGatkGenotype() {
		return gatkGenotype;
	}

}
