package dataStructure;

import algorithm.IntervalAbs;
import net.sf.samtools.SAMRecord;

public class PairedReadEntry extends FragmentReadEntry {
	static long readMateUnmappedCount = 0;
	static long readMateSplittedChromosome = 0;
	static long bothSameOrientated = 0;
	static long highDistance = 0;

	public PairedReadEntry(SAMRecord samRecord) {
		super(samRecord);
		// TODO Auto-generated constructor stub
	}

	@Override
	public SAMRecord getSamRecord() {
		// TODO Auto-generated method stub
		return super.getSamRecord();
	}

	@Override
	public void setSamRecord(SAMRecord samRecord) {
		// TODO Auto-generated method stub
		super.setSamRecord(samRecord);
	}

	@Override
	public void analyzeCoverage() {
		// in this analysis, skip and count the read pair where either read
		// or mate are flagged as unmapped
		if (this.getSamRecord().getReadUnmappedFlag()
				|| this.getSamRecord().getMateUnmappedFlag()) {
			readMateUnmappedCount++;
			readCount--;
			return;
		}

		// in this analysis, skip and count the reads located on different
		// chromosomes
		// ReferenceIndex seems to be the Index which is 0-based Chromosome
		// Number
		if (this.getSamRecord().getReferenceIndex() != this.getSamRecord()
				.getMateReferenceIndex()) {
			readMateSplittedChromosome++;
			readCount--;
			return;
		}

		// in this analysis, skip and count the pairs where read and mate
		// have same orientation
		if (this.getSamRecord().getReadNegativeStrandFlag()
				&& this.getSamRecord().getMateNegativeStrandFlag()) {
			bothSameOrientated++;
			readCount--;
			return;
		}

		if (!this.getSamRecord().getReadNegativeStrandFlag()
				&& !this.getSamRecord().getMateNegativeStrandFlag()) {
			bothSameOrientated++;
			readCount--;
			return;
		}

		// in this analysis, skip and count the read pairs with InsertSize >
		// 1000 bp
		if (Math.abs(this.getSamRecord().getInferredInsertSize()) > 1000) {
			highDistance++;
			readCount--;
			return;
		}

		// get the offset for the absolute localization
		long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
				(short) (this.getSamRecord().getReferenceIndex() + 1))
				.getOffset();
		long absAlStart = 0;
		long absAlEnd = 0;

		if (!this.getSamRecord().getReadNegativeStrandFlag()) {
			absAlStart = this.getSamRecord().getAlignmentStart() + offset;
			absAlEnd = this.getSamRecord().getAlignmentEnd() + offset;
		} else { // die Startkoordinate von reversen Reads ist im
			// BAM File nur über End-Koordinate minus Länge
			// des Reads aus dem Feld size + 1 zu berechnen
			absAlStart = offset + this.getSamRecord().getAlignmentEnd()
					- getSamRecord().getInferredInsertSize() + 1;
			absAlEnd = offset + getSamRecord().getAlignmentEnd();
		}

		// sieht so aus, als müßte ich von dem Read ein IntervalAbs
		// erzeugen, dass ich floor() übergeben kann, um das richtige
		// IntervalAbs- Objekt zurückzubekommen.

		IntervalAbs tempRead = new IntervalAbs((short) (getSamRecord()
				.getReferenceIndex() + 1), absAlEnd, absAlStart,
				(int) (absAlEnd - absAlStart), null);
		
		
	}
}
