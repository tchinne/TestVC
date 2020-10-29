function [RT_Structures, ROIBitMaps, ROI_Name] = LoadDICOMRTs(Img,Info,ImgFileName,RTSTRUCTFileName,CONTOURNumber)

RT_Structures = dicominfo([RTSTRUCTFileName],'UseVRHeuristic',false); % The strucutres are stored in the DICOM headers

% Initiate a BitMap of 64Bit
ROIBitMaps = uint64(zeros(size(Img,1),size(Img,2),size(Img,3)));

% %% Error Checking: Check whether the image corresponds to the contours delineated
% if isfield(RT_Structures.ROIContourSequence.Item_1,'ContourSequence')==0
%     disp('This RT STRUCT contains empty ROI with no contours.');
% else
%     ROI_ReferencedSOPInstanceUID = RT_Structures.ROIContourSequence.Item_1.ContourSequence.Item_1.ContourImageSequence.Item_1.ReferencedSOPInstanceUID;
%     for i = 1:numel(ImgFileName) % Go through each DICOM image and get the SOPInstanceUID
%         tempSOPInstanceUID = Info{i};
%         tempSOPInstanceUID = tempSOPInstanceUID.SOPInstanceUID;
%         SOPInstanceUID{i} = tempSOPInstanceUID;
%         
%         % Find if any SOPInstanceUID matches the ROI ReferencedSOPInstanceUID
%         NameMatch(i) = strcmp(ROI_ReferencedSOPInstanceUID,SOPInstanceUID{i});
%     end
%     
%     % Find whether NameMatch==1 (in other words there is a match
%     NameMatch = find(NameMatch==1);
%     if isempty(NameMatch) == 1 % There is no match between ROI and Image
%         warndlg('There is a DICOM Image and DICOM-RT Structure mismatch.','!!Warning!!');
%         ROIBitMask = [];
%         return
%     else
%     end
% end


%% Get PixelSpacing, ImagePositionPatient, and SliceLocation from DICOM Image Headers
for i = 1:numel(Info)
    PixelSpacing(i,:) = Info{i}.PixelSpacing;
    ImagePositionPatient(i,:) = Info{i}.ImagePositionPatient;
    if isfield(Info{i},'SliceLocation')==0 % If Slice Location does not exist
        SliceLocation(i) = Info{i}.ImagePositionPatient(3); % Pretend ImagePositionPatient is the slice location
    else
        SliceLocation(i,:) = Info{i}.SliceLocation;
    end
end

%% Get the names of the RT Structures
temp_Names = struct2cell(RT_Structures.StructureSetROISequence);

for i = 1:numel(temp_Names)
    temp = temp_Names(i);
    temp = temp{1,1}.ROIName;
    ROI_Name{i} = temp; % Name of the ROIs are now stored as ROI_Names in strings
end
NumNames = numel(ROI_Name);

%% Get the actual RT Structures
% ROIContourSequence = RT_Structures.ROIContourSequence; % Going into the actual field that contains the RT Structures
% NumRTStructures = numel(fieldnames(ROIContourSequence)); % Number of RT Structures
% Structure_ID = fieldnames(ROIContourSequence); % Get the individual RT field names

% ROIContourSequence = RT_Structures.ROIContourSequence.Item_13.ContourSequence;
% NumRTStructures = numel(fieldnames(ROIContourSequence));
% Structure_ID = fieldnames(ROIContourSequence);

%gtvIndex = find(contains(ROI_Name,'gtv','IgnoreCase',true));
% structNumber = num2cell(contourNumber);
% structureIndex = find(contains(Structure_ID,'13'));
% gtv = CONTOURName;
% ans= cell2struct(CONTOURName);
NonROIIndex = []; % For RT STRUCT with no delineated contours
%for m=1:length(gtvIndex) %NumRTStructures % Going through each ROI one at a time
    %tempStructure = getfield(ROIContourSequence,Structure_ID{i});
%     if isfield(tempStructure,'ContourSequence') == 0
%         disp('You have a RT STRUCT with no delineated contours!.')
%         NonROIIndex = [NonROIIndex;i];
%     else

%       tempStructure = tempStructure.ContourSequence;
%       NumContours = numel(fieldnames(tempStructure)); % The number of contours in each ROI
%       ContourID = fieldnames(tempStructure); % ID of the contours
        tempStructure = RT_Structures.ROIContourSequence.Item_13.ContourSequence;
        NumContours = numel(fieldnames(tempStructure)); % The number of contours in each ROI
        ContourID = fieldnames(tempStructure); % ID of the contours
        for j = 1:NumContours % Go through each contour (or slice)
            tempCoordinates = getfield(tempStructure,ContourID{j});
            tempCoordinates = tempCoordinates.ContourData; % Coordinates are in [x1;z1;y1;x2;z2;y2; etc] format
            tempCoordinates = reshape(tempCoordinates,3,numel(tempCoordinates)/3); % reshape to [x1;z1;y1,x2;z2;y2, etc] format
            
            % Determine which orientation is the slice - coordinates are the same if it is the same slice
            temp_x = range(tempCoordinates(1,:)); % First coordinate
            temp_y = range(tempCoordinates(2,:)); % Second coordinate
            temp_z = range(tempCoordinates(3,:)); % Third coordinate
            
            if temp_x == 0 % If the orientation of the contours are sagittal
                % Find out which slice we are dealing with
                temp_slice = tempCoordinates(1,:);
                tempSliceID = abs(SliceLocation-temp_slice(1));
                [~,SliceID{i}(j)] = min(tempSliceID);
                
                tempPixelSpacing = PixelSpacing(SliceID{i}(j),:);
                tempImagePositionPatient = ImagePositionPatient(SliceID{i}(j),:);
                
                temp_row = tempCoordinates(2,:)/tempPixelSpacing(1)-tempImagePositionPatient(2)/tempPixelSpacing(1)+1; % This is actually column
                temp_col = tempImagePositionPatient(3)/tempPixelSpacing(2)-tempCoordinates(3,:)/tempPixelSpacing(2)+1; % This is actually row
                
            elseif temp_y == 0 % If orientation of the contours are coronal
                % Find out which slice we are dealing with
                temp_slice = tempCoordinates(2,:);
                tempSliceID = abs(SliceLocation-temp_slice(1));
                [~,SliceID{i}(j)] = min(tempSliceID);
                
                tempPixelSpacing = PixelSpacing(SliceID{i}(j),:);
                tempImagePositionPatient = ImagePositionPatient(SliceID{i}(j),:);
                
                temp_row = tempCoordinates(1,:)/tempPixelSpacing(1)-tempImagePositionPatient(1)/tempPixelSpacing(1)+1; % This is actually column coordinates
                temp_col = tempImagePositionPatient(3)/tempPixelSpacing(2)-tempCoordinates(3,:)/tempPixelSpacing(2)+1; % This is actually row coordinates
                
            elseif temp_z ==0 % If orientation of the contours are axial
                % Find out which slice we are dealing with
                temp_slice = tempCoordinates(3,:);
                tempSliceID = abs(ImagePositionPatient(:,3)-temp_slice(1));  % Sarah changed from SliceLocation, as this code wasn't taking into acccount patient orientation in scanner
                [~,SliceID{i}(j)] = min(tempSliceID);
                
                tempPixelSpacing = PixelSpacing(SliceID{i}(j),:);
                tempImagePositionPatient = ImagePositionPatient(SliceID{i}(j),:);
                
                temp_row = tempCoordinates(1,:)/tempPixelSpacing(1)-tempImagePositionPatient(1)/tempPixelSpacing(1)+1; % This is actually column coordinates
                temp_col = tempCoordinates(2,:)/tempPixelSpacing(2)-tempImagePositionPatient(2)/tempPixelSpacing(2)+1; % This is actually row coordinates
            end
            
            %% Actually turning the ROIs into bit mask
            temp_mask = {uint64(poly2mask(temp_row,temp_col,size(Img,1),size(Img,2)))}; % Turns the contour into a mask
            temp_mask = cell2struct(temp_mask,'temp_mask',1);
            
            temp_ROIBitMaps = {ROIBitMaps(:,:,SliceID{i}(j))}; % The map to be updated
            temp_ROIBitMaps = cell2struct(temp_ROIBitMaps,'temp_ROIBitMaps',1);
            
            % This create a map of zeros in order to get the proper bit position
            temp_zeros = uint64(zeros(size(Img,1),size(Img,2)));
            temp_zeros = cell2struct({temp_zeros},'temp_zeros',1);
            temp_mask2 = arrayfun(@(x) bitset(x.temp_zeros,i,'uint64'),...
                temp_zeros,'UniformOutput',false); %Set the bit for the map to be updated
            temp_mask2 = cell2mat(temp_mask2).*cell2mat(struct2cell(temp_mask));
            temp_mask2 = cell2struct({temp_mask2},'temp_mask2',1);
            
            % Update the current map with the previous ROI bitmaps
            tempBitMask = arrayfun(@(x,y) bitor(x.temp_mask2, y.temp_ROIBitMaps),...
                temp_mask2, temp_ROIBitMaps,'UniformOutput',false);
            
            ROIBitMaps(:,:,SliceID{i}(j)) = tempBitMask{1};

        end
    %end
%end
 
%return