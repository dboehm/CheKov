package dataStructure;

import net.sf.samtools.SAMRecord;

public abstract class ReadEntry {
	private SAMRecord samRecord;
	public static long readCount = 0;

	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public void setSamRecord(SAMRecord samRecord) {
		this.samRecord = samRecord;
	}

	public abstract void analyzeCoverage();
}
