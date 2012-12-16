package dataStructure;

import net.sf.samtools.SAMRecord;

public abstract class ReadEntry {
	private static int readCount;
	private SAMRecord samRecord;
	private int softClippedBases = 0;
	private int hardClippedBases = 0;
	private int rawReadLength = 0;
	private int effReadLength = 0;
	private int deletedTaggedBases = 0;
	private int insertedTaggedBases = 0;

	
	// Getter and Setter
	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public void setSamRecord(SAMRecord samRecord) {
		this.samRecord = samRecord;
	}

	public abstract void analyzeCoverage();
	public abstract void analyseQuality();

	public int getHardClippedBases() {
		return hardClippedBases;
	}

	public void setHardClippedBases(int hardClippedBases) {
		this.hardClippedBases = hardClippedBases;
	}

	public int getSoftClippedBases() {
		return softClippedBases;
	}

	public void setSoftClippedBases(int softClippedBases) {
		this.softClippedBases = softClippedBases;
	}

	public int getRawReadLength() {
		return rawReadLength;
	}

	public void setRawReadLength(int rawReadLength) {
		this.rawReadLength = rawReadLength;
	}

	public static int getReadCount() {
		return readCount;
	}

	public static void setReadCount(int readCount) {
		ReadEntry.readCount = readCount;
	}

	public int getDeletedTaggedBases() {
		return deletedTaggedBases;
	}

	public void setDeletedTaggedBases(int deletedTaggedBases) {
		this.deletedTaggedBases = deletedTaggedBases;
	}

	public int getInsertedTaggedBases() {
		return insertedTaggedBases;
	}

	public void setInsertedTaggedBases(int insertedTaggedBases) {
		this.insertedTaggedBases = insertedTaggedBases;
	}

	public int getEffReadLength() {
		return effReadLength;
	}

	public void setEffReadLength(int effReadLength) {
		this.effReadLength = effReadLength;
	}

	public abstract void collectQualities();
}
