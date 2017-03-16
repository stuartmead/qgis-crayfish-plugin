/*
Crayfish - A collection of tools for TUFLOW and other hydraulic modelling packages
Copyright (C) 2016 Lutra Consulting

info at lutraconsulting dot co dot uk
Lutra Consulting
23 Chestnut Close
Burgess Hill
West Sussex
RH15 8HN

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
if the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "crayfish.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "crayfish_dataset.h" 
#include "crayfish_output.h"
#include "crayfish_mesh.h"

#include "crayfish_hdf5.h"
#include "contrib/tinyxml2.h"
#include "elem/crayfish_e4q.h"

//class titanSnapshot;
//class outputUpdaterTitan;
//static int readXdmfOutput(ElementOutput *o, int iT, int iF, outputUpdaterXdmf *updater);
/*
class titanSnapshot {
public:
    snapshot(){}
    double time;
    std::string connections; //String of the path to element connections
    int nConnections; //Number of connections
    int nNodes; //Number of nodes
    std::string coordinates; //String of the path to node co-ordinates
    std::vector<std::string> properties; //Name of each property
    std::vector<std::string> propertyPaths; //Path to each property (cell-centred, must be nConnections long)
};

class outputUpdaterTitan : public outputUpdater {
public:
    outputUpdaterTitan(){
        lastHdfFilename = "";
        lastDataSet = "";
        HdFile = NULL;
    }
    ~outputUpdaterTitan(){
        if (HdFile) delete HdFile;
    }
    std::vector<titanSnapshot> snaps;
    std::string lastHdfFilename, lastDataSet; //Avoid reloading if file or dataset is the same
    HdfFile *HdFile;
    //TODO here
    QVector<hsize_t> dataDim;
    QVector<double> data;
    int update(Output *o, int iTime, int iField){
        return readXdmfOutput(static_cast<ElementOutput*>(o), iTime, iField, this);
    }
};


tinyxml2::XMLElement *getCheckChild(tinyxml2::XMLElement *parent, std::string name){
    tinyxml2::XMLNode *child = parent->FirstChildElement(name.c_str());
    if (!child){
        qDebug("XML Element not found: %s", name.c_str());
        return NULL;
    }
    return child->ToElement();
}

tinyxml2::XMLElement *getCheckSibling(tinyxml2::XMLElement *from, std::string name){
    tinyxml2::XMLNode *child = from->NextSiblingElement(name.c_str());
    if (!child){
        qDebug("XML Element not found: %s", name.c_str());
        return NULL;
    }
    return child->ToElement();
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if ( item !="" )
            elems.push_back(item);
    }
    return elems;
}*/

static HdfFile openHdfFile(const QString& fileName)
{
    HdfFile file(fileName);
    if (!file.isValid())
    {
      throw LoadStatus::Err_UnknownFormat;
    }
    return file;
}

static HdfGroup openHdfGroup(const HdfFile& hdfFile, const QString& name)
{
    HdfGroup grp = hdfFile.group(name);
    if (!grp.isValid())
    {
      throw LoadStatus::Err_UnknownFormat;
    }
    return grp;
}

static HdfGroup openHdfGroup(const HdfGroup& hdfGroup, const QString& name)
{
    HdfGroup grp = hdfGroup.group(name);
    if (!grp.isValid())
    {
      throw LoadStatus::Err_UnknownFormat;
    }
    return grp;
}

static HdfDataset openHdfDataset(const HdfGroup& hdfGroup, const QString& name)
{
    HdfDataset dsFileType = hdfGroup.dataset(name);
    if (!dsFileType.isValid())
    {
      throw LoadStatus::Err_UnknownFormat;
    }
    return dsFileType;
}


static QString openHdfAttribute(const HdfFile& hdfFile, const QString& name)
{
    HdfAttribute attr = hdfFile.attribute(name);
    if (!attr.isValid())
    {
      throw LoadStatus::Err_UnknownFormat;
    }
    return attr.readString();
}

static Mesh* parseTitanMesh(HdfGroup gMesh)
{
    Mesh::Nodes nodes;
    Mesh::Elements elements;

    HdfDataset point = openHdfDataset(gMesh, "Points");

    QVector<hsize_t> nNds = point.dims();
    QVector<float> points = point.readArray();
    
    nodes.resize(nNds.at(0));
    qDebug(QString("Resized nodes to %1").arg(nNds.at(0)).toLocal8Bit().constData());
    
    for (int n = 0; n < nNds.at(0); ++n){
        nodes[n].setId(n);
        nodes[n].x = points[nNds[1]*n];
        nodes[n].y = points[nNds[1]*n+1];        
    }

    HdfDataset conns = openHdfDataset(gMesh, "Connections");

    QVector<hsize_t> nElm = conns.dims();
    QVector<int> elem_conn = conns.readArrayInt();

    elements.resize(nElm.at(0));
    qDebug(QString("Resized elements to %1").arg(nElm.at(0)).toLocal8Bit().constData());
    for (int e = 0; e < nElm.at(0); ++e)
    {
        elements[e].setId(e);
        elements[e].setEType(Element::E4Q);
        QVector<uint> idx(nElm.at(1));
        for (int fi=0; fi < nElm.at(1); ++fi)
        {
            idx[fi] = elem_conn[nElm[1]*e + fi];
        }

        elements[e].setP(idx.data());
    }

    return new Mesh(nodes, elements);

}

Mesh* Crayfish::loadTitan2D(const QString& fileName, LoadStatus* status)
{
    if (status) status->clear();
    Mesh* mesh = 0;

    try
    {
        HdfFile hdfFile = openHdfFile(fileName);
        HdfGroup geom = openHdfGroup(hdfFile, "Mesh"); //Mesh group
        HdfGroup props = openHdfGroup(hdfFile, "Properties");

        mesh = parseTitanMesh(geom);
        
        //Add pileheight
        DataSet* pHt = new DataSet(fileName);
        pHt->setType(DataSet::Scalar);
        pHt->setName("Pile Height");
        pHt->setIsTimeVarying(false);//Setting time varying to false - only read in 1 data
        
        ElementOutput* phos = new ElementOutput;
        phos->init(mesh->elements().size(), false);
        phos->time = 0.0;//Time stuff

        HdfDataset pile = openHdfDataset(props, "PILE_HEIGHT");    
        QVector<float> pileVals = pile.readArray();
        for(int i = 0; i < mesh->elements().size(); ++i)
        {
            phos->getValues()[i] = pileVals[i];     
        }
        pHt->addOutput(phos);

        //Add x, y momentum
        DataSet* pMom = new DataSet(fileName);
        pMom->setType(DataSet::Vector);
        pMom->setName("Momentum");
        pMom->setIsTimeVarying(false);

        ElementOutput* moms = new ElementOutput;
        moms->init(mesh->elements().size(), true);
        moms->time = 0.0;

        HdfDataset xmom = openHdfDataset(props, "XMOMENTUM");
        HdfDataset ymom = openHdfDataset(props, "YMOMENTUM");
        QVector<float> xVals = xmom.readArray();
        QVector<float> yVals = ymom.readArray();
        qDebug(QString("Momentum maximum size is %1 x and %2 y").arg(*std::max_element(xVals.constBegin(),xVals.constEnd()))
                                                                .arg(*std::max_element(yVals.constBegin(),yVals.constEnd()))
                                                                .toLocal8Bit().constData());

        for(int i = 0; i < mesh->elements().size(); ++i)
        {
            moms->getValuesV()[i].x = xVals[i];
            moms->getValuesV()[i].y = yVals[i];
            moms->getValues()[i] = moms->getValuesV()[i].length();
        }
        pMom->addOutput(moms);    

        //Sort out
        pHt->updateZRange();
        mesh->addDataSet(pHt);
        pMom->updateZRange();
        mesh->addDataSet(pMom);

        
    }

    catch (LoadStatus::Error error)
    {
        if (status) status->mLastError = (error);
        if (mesh) delete mesh;
        mesh = 0;
    }

    return mesh;

}

//Parse file - done
/*
static int parseTitanXdmfXml(const QString& datFileName, std::vector<titanSnapshot> &snaps){
    std::string dfName = datFileName.toUtf8().constData();
    std::string dirName = dfName.substr(0, dfName.find_last_of("/\\") + 1);
    tinyxml2::XMLDocument xmfFile;
    if (xmfFile.LoadFile(dfName.c_str()) != tinyxml2::XML_SUCCESS){
        qDebug("Cannot open xmf file");
        return -1;
    }
    tinyxml2::XMLElement *elem = xmfFile.FirstChildElement("Xdmf");
    if (!elem){
        qDebug("XML Element not found: Xdmf");
        return -1;
    }
    if (strcmp(elem->Attribute("Version"), "2.0") != 0){
        qDebug("Only Xmdf version 2.0 is supported");
        return -1;
    }
    elem = getCheckChild(elem, "Domain");
    if (!elem) return -1;
    elem = getCheckChild(elem, "Grid");
    if (!elem) return -1;
    if (strcmp(elem->Attribute("GridType"), "Collection") != 0){
        qDebug("Expecting Collection Grid Type");
        return -1;
    }
    if (strcmp(elem->Attribute("CollectionType"), "Temporal") != 0){
        qDebug("Expecting Temporal Collection Type");
        return -1;
    }
    for (tinyxml2::XMLNode* gridNod = getCheckChild(elem, "Grid"); gridNod != NULL; gridNod = gridNod->NextSiblingElement("Grid")){//Loop through each Grid
        titanSnapshot snap;
        double t;
        tinyxml2::XMLElement *timeNod = getCheckChild(gridNod->ToElement(), "Time");
        if (!timeNod) return -1;
        timeNod->QueryDoubleAttribute("Value", &t);
        snap.time = t;
        
        //Get topology
        tinyxml2::XMLElement *topology = getCheckChild(gridNod->ToElement(), "Topology");
        if (!topology) return -1;
        //Check grid is quad
        const char* attrText = nullptr;
        topology->Attribute("Type", attrText);
        std::string topoType = attrText;
        if (attrText.compare("QUADRILATERAL") != 0)
        {
            qDebug("ERROR: Topology is not quadrilateral type");
            return -1;
        }
        int elems;
        topology->QueryIntAttribute("Dimensions", &elems);
        snap.nConnections = elems;
        
        //Get connections
        tinyxml2::XMLElement *connData = getCheckChild(topology, "DataItem");
        if (!connData) return -1;
        std::string connPath = connData->GetText();
        snap.connections = connPath;

        //Get geometry (coordinates)
        tinyxml2::XMLElement *geomData = topology->NextSiblingElement("Geometry");
        tinyxml2::XMLElement *coordData = getCheckChild(geomData, "DataItem");
        if (!coordData) return -1;
        std::string dimens = coordData->Attribute("Dimensions");
        std::stringstream dimensS(dimens);
        dimensS >> snap.nNodes; //Should extract number of nodes
        std::string coordPath = coordData->GetText();
        snap.coordinates = coordPath;

        //Loop through attributes
        for (tinyxml2::XMLElement *attrElem = geomData->NextSiblingElement("Attribute"); attrElem != NULL; attrElem = attrElem->NextSiblingElement("Attribute")){
            snap.properties.push_back(attrElem->Attribute("Name"));
            tinyxml2::XMLElement *dataPropPath = getCheckChild(attrElem, "DataItem");
            if(!dataPropPath) return -1;
            std::string pPath = dataPropPath->GetText();
            snap.propertyPaths.push_back(pPath);      
        }
        snaps.push_back(snap);
    }
    return 0;
}
*/
//On the fly loading - not done.
/*
Mesh::DataSets Crayfish::loadTitanXdmfDataSet(const QString& datFileName, LoadStatus* status, Mesh* mesh) //NOT const Mesh* mesh, 
{
    outputUpdaterTitan *updater = new outputUpdaterTitan();

    cast<Mesh*>(mesh)->updater = updater;

    //get Dimensions
    if (parseTitanXdmfXml(datFileName, updater->snaps) < 0){
        if (status) status->mLastError = LoadStatus::Err_UnknownFormat;
        return 0;
    }

    std::string path = updater->snaps[0].connections;
    updater->lastHdfFilename = connections.substr(0, connections.find_last_of(":"));
    if (updater->HdFile) delete updater->HdFile;
    updater->HdFile = new HdfFile(QString::fromStdString(updater->lastHdfFilename));
    if (!updater->HdFile->isValid())
    {
        qDebug("Cannot read h5 file %s", updater->lastHdfFilename.c_str());
        if (status) status->mLastError = LoadStatus::Err_UnknownFormat;
        return 0;
    }
    //We suppose all timesteps to have the same structure
    Mesh::DataSets datasets;
    for (size_t iF = 0; iF < updater->snaps[0].names.size(); iF++){
        DataSet* ds = new DataSet(datFileName, iF);
        ds->setIsTimeVarying(true);
        ds->setName(QString::fromStdString(updater->snaps[0].names[iF]), false);
        if (updater->snaps[0].HyperSlabs[iF][2][1] == 1)
            ds->setType(DataSet::Scalar);
        else if (updater->snaps[0].HyperSlabs[iF][2][1] == 3)
            ds->setType(DataSet::Vector);
        else
        {
            qDebug("Wrong HyperSlab value: second column should terminate by 1 for scalars and 3 for vectors");
            if (status) status->mLastError = LoadStatus::Err_UnknownFormat;
            return Mesh::DataSets();
        }
        datasets.append(ds);
    }

    for (size_t iT = 0; iT < updater->snaps.size(); iT++){
        for (size_t iF = 0; iF < updater->snaps[iT].names.size(); iF++){
            ElementOutput* o = NULL;
            o = new ElementOutput(iT);
            o->time = updater->snaps[iT].time;
            datasets[iF]->addOutput(o);
        }
    }
    return datasets;
}
*/

//On the fly reading - not done.
/*
static int readXdmfOutput(ElementOutput *o, int iT, int iF, outputUpdaterXdmf *updater){
    snapshot &snap = updater->snaps[iT];
    std::string path = snap.paths[iF];
    std::string fname = path.substr(0, path.find_last_of(":"));
    if (strcmp(fname.c_str(), updater->lastHdfFilename.c_str()) != 0){
        if (updater->HdFile) delete updater->HdFile;
        updater->HdFile = new HdfFile(QString::fromStdString(fname));
        if (!updater->HdFile->isValid())
        {
            qDebug("Cannot read h5 file : %s", fname.c_str());
            return -1;
        }
        updater->lastHdfFilename = fname;
    }
    QStringList groupNames = updater->HdFile->groups();
    path = path.substr(path.find_last_of(":") + 1);
    std::vector<std::string> pathSplit = split(path, '/');
    if (!groupNames.contains(QString::fromStdString(pathSplit[0]))){
        qDebug("Datagroup not found : %s", pathSplit[0].c_str());
        return -1;
    }
    HdfGroup dataGroup = updater->HdFile->group(QString::fromStdString(pathSplit[0]));
    if (!dataGroup.isValid()){
        qDebug("Datagroup not valid");
        return -1;
    }
    for (size_t iP = 1; iP < pathSplit.size() - 1; iP++){
        if (!groupNames.contains(QString::fromStdString(pathSplit[iP]))){
            qDebug("Datagroup not found : %s", pathSplit[iP].c_str());
            return -1;
        }
        dataGroup = dataGroup.group(QString::fromStdString(pathSplit[iP]));
        if (!dataGroup.isValid()){
            qDebug("Datagroup not valid");
            return -1;
        }
    }
    QStringList dataNames = dataGroup.datasets();
    if (!dataNames.contains(QString::fromStdString(pathSplit[pathSplit.size() - 1])))
    {
        qDebug("Dataset not having required arrays : %s", pathSplit[pathSplit.size() - 1].c_str());
    }
    if (strcmp(updater->lastDataSet.c_str(), pathSplit[pathSplit.size() - 1].c_str()) != 0){
        updater->lastDataSet = pathSplit[pathSplit.size() - 1];
        HdfDataset HDFdata = dataGroup.dataset(QString::fromStdString(pathSplit[pathSplit.size() - 1]));
        updater->data = HDFdata.readArrayDouble();  //Todo should be read only once for all fields, as for file opening
        updater->dataDim = HDFdata.dims();
    }
    if (snap.HyperSlabs[iF][2][1] == 1){ //Scalar
        o->init(snap.HyperSlabs[iF][2][0], false);
        int iElem = 0;
        for (int iSlab = 0; iSlab < snap.HyperSlabs[iF][2][0]; iSlab++){
            int idxI = snap.HyperSlabs[iF][0][0] + snap.HyperSlabs[iF][1][0] * iSlab;
            int idxJ = snap.HyperSlabs[iF][0][1];
            o->getValues().data()[iElem] = updater->data[updater->dataDim[1] * idxI + idxJ];
            iElem++;
        }
    }
    else { //Vector
        o->init(snap.HyperSlabs[iF][2][0], true);
        int iElem = 0;
        ElementOutput::float2D* data = o->getValuesV().data();
        for (int iSlab = 0; iSlab < snap.HyperSlabs[iF][2][0]; iSlab++){
            for (int jSlab = 0; jSlab < snap.HyperSlabs[iF][2][1]; jSlab++){
                int idxI = snap.HyperSlabs[iF][0][0] + snap.HyperSlabs[iF][1][0] * iSlab;
                int idxJx = snap.HyperSlabs[iF][0][1];
                int idxJy = snap.HyperSlabs[iF][0][1] + snap.HyperSlabs[iF][1][1];

                data[iElem].x = updater->data[updater->dataDim[1] * idxI + idxJx];
                data[iElem].y = updater->data[updater->dataDim[1] * idxI + idxJy];
                o->getValues().data()[iElem] = data[iElem].length();
            }
            iElem++;
        }

    }
    return 0;
}
*/