package uamds.vis;

import java.awt.geom.Point2D;
import java.util.List;
import java.util.function.IntUnaryOperator;

import hageldave.imagingkit.core.Pixel;
import hageldave.jplotter.color.ColorOperations;
import hageldave.jplotter.interaction.kml.CoordSysPanning;
import hageldave.jplotter.interaction.kml.CoordSysScrollZoom;
import hageldave.jplotter.renderables.Lines;
import uamds.misc.Utils2;
import uamds.optimization.generic.numerics.MatCalc;
import uamds.other.Utils;
import uamds.other.NRV;
import uamds.other.NRVSet;
import uamds.other.Ref;

public class DistributionTravel<M> extends DistributionPlot<M> {

	MatCalc<M> mc;
	Ref<List<NRVSet<M>>> distribs = new Ref<>();
	Ref<String[]> dirLabels = new Ref<>();
	Lines meanTraces = new Lines().setStrokePattern(0b1100110011001100);
	boolean coloredLines = true;
	public IntUnaryOperator colorForDistribProvider = (i)->0xff222222;
	
	public DistributionTravel(MatCalc<M> mc) {
		super(mc,1);
		this.mc = mc;
		Utils2.addExportMenuOnMiddleClick(canvas);
		distribs.addListener(this::onDistributionChange);
		this.completeRendererBG.addItemToRender(meanTraces);
		
//		new CoordSysPanning(canvas, coordsys).register();
//		new CoordSysScrollZoom(canvas, coordsys).register();
	}
	
	public void setDistribs(List<NRVSet<M>> distribs, String[] labels) {
		this.dirLabels.set(labels);
		setDistribs(distribs);
	}
	
	public void setDistribs(List<NRVSet<M>> distribs) {
		this.distribs.set(distribs);
	}
	
	void onDistributionChange(List<NRVSet<M>> distribTraces) {
		meanTraces.removeAllSegments();
		int nDistribs = 3;
		ensureCapacity(distribTraces.size()*nDistribs);
		for(int i=0; i<distribTraces.size(); i++) {
			purgeDistribution(i*nDistribs+0);
			purgeDistribution(i*nDistribs+1);
			purgeDistribution(i*nDistribs+2);
//			purgeDistribution(i*nDistribs+3);
		}
		if(distribTraces!=null && distribTraces.size()>0) {
			// trace means
			for(int i=0; i<distribTraces.size(); i++) {
				int color = colorForDistribProvider.applyAsInt(i);
				NRVSet<M> trace = distribTraces.get(i);
				Point2D[] points = trace.stream().map(nrv->nrv.mean).map(m->new Point2D.Double(mc.get(m, 0),mc.get(m, 1))).toArray(Point2D[]::new);
				meanTraces.addLineStrip(points).forEach(seg->seg.setColor(color));
			}
			for(int i=0; i<distribTraces.size(); i++) {
				int color = colorForDistribProvider.applyAsInt(i);
				NRVSet<M> trace = distribTraces.get(i);
				NRV<M> t_0 = trace.get(0)			.scaleCov(0.01);
				NRV<M> t_1 = trace.get(trace.size()/2).scaleCov(0.01);
//				NRV<M> t_2 = trace.get(trace.size()*2/3).scaleCov(0.01);
				NRV<M> t_3 = trace.get(trace.size()-1).scaleCov(0.01);
				updateDistribution(i*nDistribs+0, t_0, ColorOperations.scaleColorAlpha(color, 0.4), "");
				updateDistribution(i*nDistribs+1, t_1, ColorOperations.scaleColorAlpha(color, 0.7), "");
//				updateDistribution(i*nDistribs+2, t_2, ColorOperations.scaleColorAlpha(color, 0.7), "");
				updateDistribution(i*nDistribs+2, t_3, color, "");
			}
//			if(!dirLabels.isNull()) {
//				for(int i=0; i<distribTraces.get(0).length; i++) {
//					Text txt = new Text(dirLabels.r[i], 12, Font.PLAIN);
//					txt.setOrigin(new Point2D.Double(distribTraces.get(0)[i][0], distribTraces.get(0)[i][1]));
//					forgrnd.addItemToRender(txt);
//				}
//			}
		}
		canvas.scheduleRepaint();
	}
	
}
