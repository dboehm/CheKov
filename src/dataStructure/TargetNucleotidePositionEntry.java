package dataStructure;

import reference.ReferenceReadPosition;
import datatypes.AlterationType;

/* this is an entry for a single bp targeted genome position out of the target, 
 * which shows any changes of any read on this position
 */

public class TargetNucleotidePositionEntry implements
		Comparable<TargetNucleotidePositionEntry> {
	// quite another job to obtain the coverage at all, and the
	// refAllelReadCount
	private Integer coverage = 1;
	private Integer refAllelReadCount = 0;

	private Integer altAllelReadCount = 1;
	private char refAllel;
	private char altAllel;
	private AlterationType alterationTypeAtPosition;
	private Long posInAbsGenome;
	private Integer homoPolymerLength;
	private ReferenceReadPosition referenceReadPostion;

	// Constructor
	public TargetNucleotidePositionEntry(long posInAbsGenome,
			AlterationType alterationTypeAtPosition, char refAllel,
			char altAllel, ReferenceReadPosition referenceReadPostion) {
		this.posInAbsGenome = posInAbsGenome;
		this.alterationTypeAtPosition = alterationTypeAtPosition;
		this.refAllel = refAllel;
		this.altAllel = altAllel;
		this.referenceReadPostion = referenceReadPostion;
	}

	public Integer getCoverage() {
		return coverage;
	}

	public void setCoverage(Integer coverage) {
		this.coverage = coverage;
	}

	public Integer getAltAllelReadCount() {
		return altAllelReadCount;
	}

	public void setAltAllelReadCount(Integer altAllelReadCount) {
		this.altAllelReadCount = altAllelReadCount;
	}

	public AlterationType getAlterationTypeAtPosition() {
		return alterationTypeAtPosition;
	}

	public void setAlterationTypeAtPosition(
			AlterationType alterationTypeAtPosition) {
		this.alterationTypeAtPosition = alterationTypeAtPosition;
	}

	public Integer getHomoPolymerLength() {
		return homoPolymerLength;
	}

	public void setHomoPolymerLength(Integer homoPolymerLength) {
		this.homoPolymerLength = homoPolymerLength;
	}

	public Long getPosInAbsGenome() {
		return posInAbsGenome;
	}

	public void setPosInAbsGenome(Long posInAbsGenome) {
		this.posInAbsGenome = posInAbsGenome;
	}

	public Integer getRefAllelReadCount() {
		return refAllelReadCount;
	}

	public void setRefAllelReadCount(Integer refAllelReadCount) {
		this.refAllelReadCount = refAllelReadCount;
	}

	public char getRefAllel() {
		return refAllel;
	}

	public void setRefAllel(char refAllel) {
		this.refAllel = refAllel;
	}

	public char getAltAllel() {
		return altAllel;
	}

	public void setAltAllel(char altAllel) {
		this.altAllel = altAllel;
	}

	public ReferenceReadPosition getReferenceReadPostion() {
		return referenceReadPostion;
	}

	public void setReferenceReadPostion(
			ReferenceReadPosition referenceReadPostion) {
		this.referenceReadPostion = referenceReadPostion;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof TargetNucleotidePositionEntry) {
			TargetNucleotidePositionEntry temp = (TargetNucleotidePositionEntry) obj;
			return (temp.posInAbsGenome == this.posInAbsGenome);
		}
		return false;

	}

	@Override
	public int compareTo(TargetNucleotidePositionEntry o) {
		if (this.equals(o))
			return 0;
		/*
		 * the range is from ~3E9 which is to much for an int, /10 maximizes
		 * values
		 */
		return (int) (this.posInAbsGenome / 10 - o.posInAbsGenome / 10);
	}

	public int calculateHomoPolymerLength() {

		return altAllel;

	}

}
