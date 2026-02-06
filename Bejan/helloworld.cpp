#include <iostream>
#include <vector>
#include <numeric>   
#include <cmath>    
#include <algorithm> 
#include <fstream>   
#include <map>       
#include <limits>    
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"
#include <DGtal/geometry/curves/ArithmeticalDSSComputer.h>
#include <DGtal/geometry/curves/GreedySegmentation.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace DGtal;
using namespace Z2i;

struct GrainMetrics {
    int id;
    double area;
    double perimeter;
    double circularity;
};

struct Stats {
    double mean = 0.0, stdev = 0.0, min = 0.0, max = 0.0;
};

struct RiceProfile {
    std::string name;
    double meanCirc;
    double stdCirc;
    Color color; 
};

template<class T>
Curve getBoundary(T & object)
{
    const DigitalSet& set = object.pointSet();
    if (set.empty()) return Z2i::Curve();

    Point lower = set.domain().upperBound(); 
    Point upper = set.domain().lowerBound();

    for (const auto& p : set) {
        lower = lower.inf(p);
        upper = upper.sup(p);
    }

    KSpace kSpace;
    kSpace.init(lower - Point(2,2), upper + Point(2,2), true);

    Point pIn  = *set.begin();
    Point pOut = kSpace.lowerBound(); 
    SCell s = Surfaces<KSpace>::findABel(kSpace, set, pIn, pOut);

    std::vector<SCell> boundarySurfels;
    SurfelAdjacency<2> sAdj(true);
    Surfaces<KSpace>::track2DBoundary(boundarySurfels, kSpace, sAdj, set, s);

    Z2i::Curve boundaryCurve;
    boundaryCurve.initFromSCellsVector(boundarySurfels);

    return boundaryCurve;
}

Stats computeStats(const std::vector<double>& data) {
    Stats s;
    if (data.empty()) return s;
    
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    s.mean = sum / data.size();
    
    double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
    s.stdev = std::sqrt(sq_sum / data.size() - s.mean * s.mean);
    
    s.min = *std::min_element(data.begin(), data.end());
    s.max = *std::max_element(data.begin(), data.end());
    return s;
}

double getZScore(double circularity, const RiceProfile& p) {
    if (p.stdCirc == 0) return std::numeric_limits<double>::max();
    return std::abs(circularity - p.meanCirc) / p.stdCirc;
}

int main(int argc, char** argv)
{
    setlocale(LC_NUMERIC, "us_US"); 
    
    typedef ImageSelector<Domain, unsigned char>::Type Image; 
    typedef DigitalSetSelector<Domain, BIG_DS | HIGH_BEL_DS>::Type DigitalSet; 
    typedef Object<DT4_8, DigitalSet> ObjectType; 
    typedef PointVector<2,int> Point;
    typedef Z2i::Curve::PointsRange::ConstCirculator MyCirculator;
    typedef ArithmeticalDSSComputer<MyCirculator, int, 4> DSS4;
    typedef GreedySegmentation<DSS4> Decomposition4;

    std::string pathPrefix = "../RiceGrains/";

    std::vector<std::string> filenames;
    filenames.push_back("Rice_japonais_seg_bin.pgm");
    filenames.push_back("Rice_basmati_seg_bin.pgm");
    filenames.push_back("Rice_camargue_seg_bin.pgm");
    
    std::vector<std::string> riceTypes;
    riceTypes.push_back("Japonais");
    riceTypes.push_back("Basmati");
    riceTypes.push_back("Camargue");

    std::ofstream csvSummary("rice_grains_analysis.csv");
    csvSummary << "GrainType,AreaMean,AreaStdDev,AreaMin,AreaMax,PeriMean,PeriStdDev,PeriMin,PeriMax,CircMean,CircStdDev,CircMin,CircMax\n";

    std::ofstream csvDetails("rice_grains_details.csv");
    csvDetails << "GrainType,ID,Area,Perimeter,Circularity,IsOutlier\n";

    for(size_t i = 0; i < filenames.size(); i++) {
        std::string fullPath = pathPrefix + filenames[i];
        std::string grainType = riceTypes[i];
        
        std::cout << "\nAnalyzing: " << grainType << " (" << fullPath << ") ---" << std::endl;
        
        try {
            Image image = PGMReader<Image>::importPGM(fullPath); 
            DigitalSet set2d(image.domain());
            SetFromImage<DigitalSet>::append<Image>(set2d, image, 1, 255);

            std::vector<ObjectType> objects; 
            std::back_insert_iterator<std::vector<ObjectType>> inserter(objects); 
            ObjectType set2d_obj(dt4_8, set2d);
            set2d_obj.writeComponents(inserter);

            // Filter boundary grains
            std::vector<ObjectType> completeObjects;
            Domain domain = image.domain();
            Point pMin = domain.lowerBound();
            Point pMax = domain.upperBound();

            for (const auto& obj : objects) {
                bool touches = false;
                for (auto const& p : obj.pointSet()) {
                    if (p[0] == pMin[0] || p[0] == pMax[0] || p[1] == pMin[1] || p[1] == pMax[1]) {
                        touches = true; break;
                    }
                }
                if (!touches) completeObjects.push_back(obj);
            }

            std::cout << "  > Grains found: " << completeObjects.size() << endl;

            std::vector<GrainMetrics> rawMetrics;

            for (size_t k = 0; k < completeObjects.size(); k++) {
                try {
                    Z2i::Curve boundaryCurve = getBoundary(completeObjects[k]);
                    Z2i::Curve::PointsRange range = boundaryCurve.getPointsRange();
                    Decomposition4 theDecomposition(range.c(), range.c(), DSS4());
                    
                    std::vector<Point> vertices;
                    for (auto it = theDecomposition.begin(); it != theDecomposition.end(); ++it) {
                        if (vertices.empty() || vertices.back() != *it->begin()) {
                            vertices.push_back(*it->begin());
                        }
                    }
                    
                    double polygonalArea = 0.0;
                    int n = vertices.size();
                    for (int j = 0; j < n; j++) {
                        int idx = (j + 1) % n;
                        polygonalArea += (double)vertices[j][0] * vertices[idx][1];
                        polygonalArea -= (double)vertices[idx][0] * vertices[j][1];
                    }
                    polygonalArea = std::abs(polygonalArea) / 2.0;
                    
                    double polygonalPerimeter = 0.0;
                    for (int j = 0; j < n; j++) {
                        int idx = (j + 1) % n;
                        double dx = vertices[j][0] - vertices[idx][0];
                        double dy = vertices[j][1] - vertices[idx][1];
                        polygonalPerimeter += std::sqrt(dx*dx + dy*dy);
                    }

                    double circularity = 0.0;
                    if (polygonalPerimeter > 0.0) {
                        circularity = (4.0 * M_PI * polygonalArea) / (polygonalPerimeter * polygonalPerimeter);
                    }

                    rawMetrics.push_back({(int)k, polygonalArea, polygonalPerimeter, circularity});

                } catch (const DGtal::InputException& e) { continue; }
            }

            // Stats & CSV
            std::vector<double> areasForSort;
            for (const auto& m : rawMetrics) areasForSort.push_back(m.area);
            std::sort(areasForSort.begin(), areasForSort.end());

            double lowerBound = -1.0, upperBound = 999999.0;
            if (!areasForSort.empty()) {
                double q1 = areasForSort[areasForSort.size() / 4];
                double q3 = areasForSort[areasForSort.size() * 3 / 4];
                double iqr = q3 - q1;
                lowerBound = q1 - 1.5 * iqr;
                upperBound = q3 + 1.5 * iqr;
            }

            std::vector<double> cleanAreas, cleanPeris, cleanCircs;
            for (const auto& m : rawMetrics) {
                bool isOutlier = (m.area < lowerBound || m.area > upperBound);
                csvDetails << grainType << "," << m.id << "," << m.area << "," << m.perimeter << "," << m.circularity << "," << (isOutlier ? 1 : 0) << "\n";

                if (!isOutlier) {
                    cleanAreas.push_back(m.area);
                    cleanPeris.push_back(m.perimeter);
                    cleanCircs.push_back(m.circularity);
                }
            }

            Stats sArea = computeStats(cleanAreas);
            Stats sPeri = computeStats(cleanPeris);
            Stats sCirc = computeStats(cleanCircs);

            std::cout << "  > Avg Area: " << sArea.mean << " | Avg Circ: " << sCirc.mean << std::endl;

            csvSummary << grainType << ","
                       << sArea.mean << "," << sArea.stdev << "," << sArea.min << "," << sArea.max << ","
                       << sPeri.mean << "," << sPeri.stdev << "," << sPeri.min << "," << sPeri.max << ","
                       << sCirc.mean << "," << sCirc.stdev << "," << sCirc.min << "," << sCirc.max << "\n";

            // Visualization (Single Grain Analysis)
            if (!completeObjects.empty()) {
                Z2i::Curve c = getBoundary(completeObjects[0]);
                Z2i::Curve::PointsRange range = c.getPointsRange();
                Decomposition4 theDecomposition( range.c(), range.c(), DSS4() );

                Board2D aBoard;
                aBoard << completeObjects[0]; 
                for (auto it = theDecomposition.begin(); it != theDecomposition.end(); ++it){
                    aBoard << SetMode( "ArithmeticalDSS", "Points" )
                           << it->primitive()
                           << SetMode( "ArithmeticalDSS", "BoundingBox" )
                           << CustomStyle("ArithmeticalDSS/BoundingBox", new CustomPenColor(Color::Red))
                           << it->primitive();
                }
                std::string pdfName = "Analyzed_" + grainType + ".pdf";
                aBoard.saveCairo(pdfName.c_str(), Board2D::CairoPDF);
            }
        } catch (...) { std::cerr << "Failed to read " << fullPath << std::endl; }
    }
    csvSummary.close();
    csvDetails.close();

    std::vector<RiceProfile> profiles;
    profiles.push_back({"Japonais", 0.858839, 0.0210094, Color::Red});
    profiles.push_back({"Basmati",  0.525403, 0.0410983, Color::Green});
    profiles.push_back({"Camargue", 0.727868, 0.0441186, Color::Blue});

    std::vector<std::string> mixedFiles;
    mixedFiles.push_back("Rice_mixed2_seg_bin.pgm");
    mixedFiles.push_back("Rice_mixed3_seg_bin.pgm");

    for (const auto& file : mixedFiles) {
        std::string fullPath = pathPrefix + file;
        std::cout << "\nClassifying: " << file << " ---" << std::endl;

        try {
            Image image = PGMReader<Image>::importPGM(fullPath); 
            DigitalSet set2d(image.domain());
            SetFromImage<DigitalSet>::append<Image>(set2d, image, 1, 255);

            std::vector<ObjectType> objects; 
            std::back_insert_iterator<std::vector<ObjectType>> inserter(objects); 
            ObjectType set2d_obj(dt4_8, set2d);
            set2d_obj.writeComponents(inserter);

            Board2D board;
            std::map<std::string, int> counts;

            Domain domain = image.domain();
            Point pMin = domain.lowerBound(), pMax = domain.upperBound();

            for (const auto& obj : objects) {
                bool touches = false;
                for (auto const& p : obj.pointSet()) {
                    if (p[0] == pMin[0] || p[0] == pMax[0] || p[1] == pMin[1] || p[1] == pMax[1]) {
                        touches = true; break;
                    }
                }
                if (touches) continue;

                Z2i::Curve boundaryCurve = getBoundary(obj);
                Z2i::Curve::PointsRange range = boundaryCurve.getPointsRange();
                Decomposition4 theDecomposition(range.c(), range.c(), DSS4());
                
                std::vector<Point> vertices;
                for (auto it = theDecomposition.begin(); it != theDecomposition.end(); ++it) {
                    if (vertices.empty() || vertices.back() != *it->begin()) vertices.push_back(*it->begin());
                }
                
                double area = 0.0, perimeter = 0.0;
                int n = vertices.size();
                for (int j = 0; j < n; j++) {
                    int idx = (j + 1) % n;
                    area += (double)vertices[j][0] * vertices[idx][1] - (double)vertices[idx][0] * vertices[j][1];
                    double dx = vertices[j][0] - vertices[idx][0];
                    double dy = vertices[j][1] - vertices[idx][1];
                    perimeter += std::sqrt(dx*dx + dy*dy);
                }
                area = std::abs(area) / 2.0;
                double circularity = (perimeter > 0) ? (4.0 * M_PI * area) / (perimeter * perimeter) : 0;

                double bestScore = std::numeric_limits<double>::max();
                int bestProfileIndex = -1;

                for(size_t i=0; i<profiles.size(); ++i) {
                    double zScore = getZScore(circularity, profiles[i]);
                    if (zScore < bestScore) {
                        bestScore = zScore;
                        bestProfileIndex = i;
                    }
                }

                if (bestProfileIndex != -1) {
                    const RiceProfile& best = profiles[bestProfileIndex];
                    counts[best.name]++;
                    
                    board.setPenColor(best.color); 
                    for (auto const& p : obj.pointSet()) {
                        board << p;
                    }
                }
            }

            std::cout << "  > Classification Counts:" << std::endl;
            for(auto const& [name, count] : counts) {
                std::cout << "    - " << name << ": " << count << std::endl;
            }

            std::string outName = "Classified_" + file + ".pdf";
            board.saveCairo(outName.c_str(), Board2D::CairoPDF);
            std::cout << "  > Visualization saved: " << outName << std::endl;

        } catch (...) { std::cerr << "Failed to read " << fullPath << std::endl; }
    }

    std::cout << "\nTasks complete." << std::endl;
    return 0;
}