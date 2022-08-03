package uamds.misc;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.MenuItem;
import java.awt.PopupMenu;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;
import java.util.function.DoubleBinaryOperator;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.apache.pdfbox.pdmodel.PDDocument;

import hageldave.imagingkit.core.Img;
import hageldave.imagingkit.core.io.ImageSaver;
import hageldave.imagingkit.core.io.ImageSaver.ImageSaverException;
import hageldave.jplotter.canvas.JPlotterCanvas;
import hageldave.jplotter.color.ColorMap;
import hageldave.jplotter.interaction.kml.CoordSysViewSelector;
import hageldave.jplotter.pdf.PDFUtils;
import hageldave.jplotter.renderers.CoordSysRenderer;
import hageldave.jplotter.svg.SVGUtils;
import uamds.optimization.generic.numerics.MatCalc;

public class Utils2 {

	public static double[][][] sampling2DCoordinates(double xmin, double xmax, double ymin, double ymax, int resolution){
		double[][] X = new double[resolution][resolution];
		double[][] Y = new double[resolution][resolution];
		for(int j=0; j<X.length;j++) {
			for(int i=0; i<X[0].length; i++) {
				double x = xmin + i*(xmax-xmin)/(resolution-1);
				double y = ymin + j*(ymax-ymin)/(resolution-1);
				X[j][i] = x;
				Y[j][i] = y;
			}
		}
		return new double[][][] {X,Y};
	}
	
	public static double[][] sample2D(DoubleBinaryOperator fn, double[][] X, double[][] Y){
		double[][] Z = new double[X.length][X[0].length];
		for(int j=0; j<X.length;j++) {
			for(int i=0; i<X[0].length; i++) {
				double x = X[j][i];
				double y = Y[j][i];
				Z[j][i] = fn.applyAsDouble(x, y);
			}
		}
		return Z;
	}
	
	public static double signedRoot(double order, double value) {
		double sign = Math.signum(value);
		return sign * Math.pow(value*sign, 1.0/order);
	}
	
	public static <M> Img mat2img(MatCalc<M> mc, M mat, double min, double max, ColorMap cmap) {
		Img img = new Img(mc.numCols(mat), mc.numRows(mat));
		double divByRange = 1.0/(max-min);
		img.forEach(px->{
			double v = mc.get(mat, px.getY(), px.getX());
			v = (v-min)*divByRange;
			v = hageldave.jplotter.util.Utils.clamp(0, v, 1);
			px.setValue(cmap.interpolate(v));
		});
		return img;
	}
	
	public static String timeNowString() {
		return LocalDateTime.now().truncatedTo(ChronoUnit.SECONDS).toString().replace(':', '-');
	}
	
	public static void exportSVG(JPlotterCanvas canvas, String filename) {
		var svg = SVGUtils.containerToSVG(canvas.asComponent().getParent());
		SVGUtils.documentToXMLFile(svg, new File(filename + ".svg"));
	}
	
	public static void exportPNG(JPlotterCanvas canvas, String filename) {
		Container container = canvas.asComponent().getParent();
		Img img = new Img(container.getSize());
		img.paint(container::paintAll);
//		Img img = canvas.toImg();
		ImageSaver.saveImage(img.toBufferedImage(),filename +".png");
	}
	
	public static void exportPDF(JPlotterCanvas canvas, String filename) {
		try {
			PDDocument doc = PDFUtils.containerToPDF(canvas.asComponent().getParent());
			doc.save(filename +".pdf");
			doc.close();
			System.out.println("exported pdf");
		} catch (IOException ex) {
			ex.printStackTrace();
		}
	}
	
	public static void exportAllFormats(JPlotterCanvas canvas, String filename) {
		exportPNG(canvas, filename);
		exportSVG(canvas, filename);
		exportPDF(canvas, filename);
	}
	
	public static void addExportMenuOnMiddleClick(JPlotterCanvas canvas) {
		PopupMenu menu = new PopupMenu();
		canvas.asComponent().add(menu);
		MenuItem svgExport = new MenuItem("SVG export");
		menu.add(svgExport);
		svgExport.addActionListener(e->{
			var svg = SVGUtils.containerToSVG(canvas.asComponent().getParent());
			SVGUtils.documentToXMLFile(svg, new File("export_"+ timeNowString() +".svg"));
			System.out.println("exported svg");
		});
		MenuItem pdfExport = new MenuItem("PDF export");
		menu.add(pdfExport);
		pdfExport.addActionListener(e->{
			try {
				PDDocument doc = PDFUtils.containerToPDF(canvas.asComponent().getParent());
				doc.save("export_"+ timeNowString() +".pdf");
				doc.close();
				System.out.println("exported pdf");
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		});
		MenuItem pngExport = new MenuItem("PNG export");
		menu.add(pngExport);
		pngExport.addActionListener(e->{
			try {
				Container container = canvas.asComponent().getParent();
				Img img = new Img(container.getSize());
				img.paint(container::paintAll);
//				Img img = canvas.toImg();
				ImageSaver.saveImage(img.toBufferedImage(), "export_"+ timeNowString() +".png");
				System.out.println("exported png");
			} catch (ImageSaverException ex) {
				ex.printStackTrace();
			}
		});
		
		canvas.asComponent().addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				if(SwingUtilities.isMiddleMouseButton(e))
					menu.show(canvas.asComponent(), e.getX(), e.getY());
			}
		});
	}
	
	public static void installFocusRequestByMouse(Component comp) {
		MouseAdapter mouseAdapter = new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				requestFocus();
			}
			@Override
			public void mousePressed(MouseEvent e) {
				requestFocus();
			}
			@Override
			public void mouseEntered(MouseEvent e) {
				requestFocus();
			}
			
			void requestFocus() {
				if(!comp.hasFocus())
					comp.requestFocus();
			}
		};
		comp.addMouseListener(mouseAdapter);
		comp.addMouseMotionListener(mouseAdapter);
	}
	
	public static CoordSysViewSelector createAreaSelectionZoom(JPlotterCanvas canvas, CoordSysRenderer coordsys) {
		return new CoordSysViewSelector(canvas, coordsys){
			@Override
			public void areaSelected(double minx, double miny, double maxx, double maxy) {
				coordsys.setCoordinateView(minx, miny, maxx, maxy);
			}
		};
	}
	
	public static JFrame createJFrameWithBoilerPlate(String title) {
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(new BorderLayout());
		frame.setTitle(title);
		frame.setPreferredSize(new Dimension(400, 400));
		return frame;
	}
	
}
