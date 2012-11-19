package dataStructure;

import algorithm.CheKov;
import algorithm.IntervalAbs;
import net.sf.samtools.SAMRecord;

public class PairedReadEntry extends FragmentReadEntry {
	private static long readMateUnmappedCount = 0;
	private static long readMateSplittedChromosome = 0;
	private static long bothSameOrientated = 0;
	private static long highDistance = 0;
	private static int pairedEndReadCount;

	public PairedReadEntry(SAMRecord samRecord) {
		super(samRecord);
		// TODO Auto-generated constructor stub
	}

	// Getter and Setter
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

	public static long getReadMateUnmappedCount() {
		return readMateUnmappedCount;
	}

	public static void setReadMateUnmappedCount(long readMateUnmappedCount) {
		PairedReadEntry.readMateUnmappedCount = readMateUnmappedCount;
	}

	public static long getReadMateSplittedChromosome() {
		return readMateSplittedChromosome;
	}

	public static void setReadMateSplittedChromosome(long readMateSplittedChromosome) {
		PairedReadEntry.readMateSplittedChromosome = readMateSplittedChromosome;
	}

	public static long getBothSameOrientated() {
		return bothSameOrientated;
	}

	public static void setBothSameOrientated(long bothSameOrientated) {
		PairedReadEntry.bothSameOrientated = bothSameOrientated;
	}

	public static long getHighDistance() {
		return highDistance;
	}

	public static void setHighDistance(long highDistance) {
		PairedReadEntry.highDistance = highDistance;
	}
	
	public static int getPairedEndReadCount() {
		return pairedEndReadCount;
	}

	public static void setPairedEndReadCount(int pairedEndReadCount) {
		PairedReadEntry.pairedEndReadCount = pairedEndReadCount;
	}

	@Override
	public void analyzeCoverage() {
		super.analyzeCoverage();
		// in this analysis, skip and count the read pair where either read
		// or mate are flagged as unmapped
		if (this.getSamRecord().getReadUnmappedFlag()
				|| this.getSamRecord().getMateUnmappedFlag()) {
			readMateUnmappedCount++;
			return;
		}

		// in this analysis, skip and count the reads located on different
		// chromosomes
		// ReferenceIndex seems to be the Index which is 0-based Chromosome
		// Number
		if (this.getSamRecord().getReferenceIndex() != this.getSamRecord()
				.getMateReferenceIndex()) {
			readMateSplittedChromosome++;
			return;
		}

		// in this analysis, skip and count the pairs where read and mate
		// have same orientation
		if (this.getSamRecord().getReadNegativeStrandFlag()
				&& this.getSamRecord().getMateNegativeStrandFlag()) {
			bothSameOrientated++;
			return;
		}

		if (!this.getSamRecord().getReadNegativeStrandFlag()
				&& !this.getSamRecord().getMateNegativeStrandFlag()) {
			bothSameOrientated++;
			return;
		}

		// in this analysis, skip and count the read pairs with InsertSize >
		// 1000 bp
		if (Math.abs(this.getSamRecord().getInferredInsertSize()) > 1000) {
			highDistance++;
			return;
		}

		// get the offset for the absolute localization
		long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
				(short) (this.getSamRecord().getReferenceIndex() + 1))
				.getOffset();
		// System.out.println("OFFSET " + offset);
		long absAlStart = 0;
		long absAlEnd = 0;

		if (!this.getSamRecord().getReadNegativeStrandFlag()) {
			absAlStart = this.getSamRecord().getAlignmentStart() + offset;
			absAlEnd = this.getSamRecord().getAlignmentEnd() + offset;
		} else { // auch bei reversen Reads ist absAlStart < absAlEnd
			absAlEnd = offset + this.getSamRecord().getAlignmentEnd();
			absAlStart = offset + getSamRecord().getAlignmentStart();
		}
		// System.out.println("READ " + getSamRecord().getReadName() +
		// " IS NEG "
		// + this.getSamRecord().getReadNegativeStrandFlag() + " ALSTART "
		// + absAlStart + " ALEND " + absAlEnd);

		// sieht so aus, als müßte ich von dem Read ein IntervalAbs
		// erzeugen, dass ich floor() übergeben kann, um das richtige
		// IntervalAbs- Objekt zurückzubekommen.
		IntervalAbs tempRead = new IntervalAbs((short) (getSamRecord()
				.getReferenceIndex() + 1), absAlEnd, absAlStart,
				(int) (absAlEnd - absAlStart), null);

		IntervalAbs floorInterval = CheKov.getIntervalTreeSet().floor(tempRead);
		if (floorInterval != null) {
			if (floorInterval.readHitsAbsInterval(tempRead.getEndAbs(),
					tempRead.getStartAbs())) {
				// System.out.println(floorInterval.getCoverage());
				return;
			}
		}
	} // end analyzeCoverage()

	@Override
	public void analyseQuality() {
		// TODO Auto-generated method stub
		super.analyseQuality();
	}

	
	
	
} // end class
