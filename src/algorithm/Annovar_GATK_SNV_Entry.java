package algorithm;

import dataStructure.GATK_SNV_entry;

public class Annovar_GATK_SNV_Entry extends GATK_SNV_entry {
	private String annotationType;
	private String gene;
	private GATK_SNV_entry gatk_SNV_entry;

	public Annovar_GATK_SNV_Entry(String annotationType, String gene,
			GATK_SNV_entry gatk_SNV_entry) {
		// this seems to be redundant, but at the moment necessary because of
		// the missing super() constructor. Needs to be addressed
		super(gatk_SNV_entry.getChr(), gatk_SNV_entry.getGenomPosition(),
				gatk_SNV_entry.getAbsGenomPosition(), gatk_SNV_entry
						.getChr_nr());
		this.setAnnotationType(annotationType);
		this.setGene(gene);
		this.setGatk_SNV_entry(gatk_SNV_entry);
	}

	public String getAnnotationType() {
		return annotationType;
	}

	public void setAnnotationType(String annotationType) {
		this.annotationType = annotationType;
	}

	public String getGene() {
		return gene;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public GATK_SNV_entry getGatk_SNV_entry() {
		return gatk_SNV_entry;
	}

	public void setGatk_SNV_entry(GATK_SNV_entry gatk_SNV_entry) {
		this.gatk_SNV_entry = gatk_SNV_entry;
	}

}
