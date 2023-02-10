clear all
load("intrinsic.mat")

filename = 'intrinsic.xml';

IntrinsicMatric = cameraParams4.IntrinsicMatrix;
Distortion = [cameraParams4.RadialDistortion, cameraParams4.TangentialDistortion];

docNode = com.mathworks.xml.XMLUtils.createDocument('Camera_Intrinsic');
docRootNode = docNode.getDocumentElement;

%Write time
dataNode = docNode.createElement('WRITE_TIME');
dataNode.appendChild(docNode.createTextNode(datestr(now)));
docRootNode.appendChild(dataNode);


% IntrinsicMatrix
dataNode = docNode.createElement('IntrinsicMatrix');
dataNode.setAttribute('type_id','opencv-matrix')

rows=docNode.createElement('rows');
rows.appendChild(docNode.createTextNode(sprintf('%d',3)));
dataNode.appendChild(rows)

cols=docNode.createElement('cols');
cols.appendChild(docNode.createTextNode(sprintf('%d',3)));
dataNode.appendChild(cols)

dt=docNode.createElement('dt');
dt.appendChild(docNode.createTextNode('d'));
dataNode.appendChild(dt)

data=docNode.createElement('data');
data.appendChild(docNode.createTextNode(sprintf('%f ',IntrinsicMatric)));
dataNode.appendChild(data)
docRootNode.appendChild(dataNode)

% DistortionCoef
dataNode = docNode.createElement('DistortionCoef');
dataNode.setAttribute('type_id','opencv-matrix')

rows=docNode.createElement('rows');
rows.appendChild(docNode.createTextNode(sprintf('%d',1)));
dataNode.appendChild(rows)

cols=docNode.createElement('cols');
cols.appendChild(docNode.createTextNode(sprintf('%d',5)));
dataNode.appendChild(cols)

dt=docNode.createElement('dt');
dt.appendChild(docNode.createTextNode('d'));
dataNode.appendChild(dt)

data=docNode.createElement('data');
data.appendChild(docNode.createTextNode(sprintf('%f ',Distortion)));
dataNode.appendChild(data)
docRootNode.appendChild(dataNode)


xmlwrite(filename, docNode)



