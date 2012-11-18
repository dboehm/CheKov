package dataStructure;

import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;

import net.sf.samtools.SAMFileHeader;

public class ProtonFileHeader extends SAMFileHeader{
SAMFileHeader sfh;
	public ProtonFileHeader(SAMFileHeader sfh) {
		this.sfh = sfh;
		Set<Entry<String, String>> attributes = sfh.getAttributes();
		System.out.println(attributes);
		System.out.println(sfh.getAttribute(CURRENT_VERSION));
		System.out.println(sfh.getAttribute(GROUP_ORDER_TAG));
		System.out.println(sfh.getAttribute(SORT_ORDER_TAG));
		System.out.println(sfh.getAttribute(VERSION_TAG));
		System.out.println(sfh.getCreator());
		
		
	}
	

}
