function varargout = computeOrbitGUI(varargin)
% COMPUTEORBITGUI MATLAB code for computeOrbitGUI.fig
%      COMPUTEORBITGUI, by itself, creates a new COMPUTEORBITGUI or raises the existing
%      singleton*.
%
%      H = COMPUTEORBITGUI returns the handle to a new COMPUTEORBITGUI or the handle to
%      the existing singleton*.
%
%      COMPUTEORBITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPUTEORBITGUI.M with the given input arguments.
%
%      COMPUTEORBITGUI('Property','Value',...) creates a new COMPUTEORBITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before computeOrbitGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to computeOrbitGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help computeOrbitGUI

% Last Modified by GUIDE v2.5 05-Aug-2017 13:59:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @computeOrbitGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @computeOrbitGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before computeOrbitGUI is made visible.
function computeOrbitGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to computeOrbitGUI (see VARARGIN)

% Choose default command line output for computeOrbitGUI
handles.output = hObject;
    
% Update handles structure
guidata(hObject, handles);

if ~isempty(varargin)
    aux = [0,varargin{:}];
    set(handles.degree,'String',num2str(aux{2}));
    set(handles.order,'String',num2str(aux{3}));
    set(handles.sun,'Value',aux{4});
    set(handles.earth,'Value',aux{5});
    set(handles.starttime,'String',aux{6});
    set(handles.frame,'String',aux{7});
    set(handles.ctype,'String',aux{8});
    set(handles.value1,'String',num2str(aux{9}));
    set(handles.value2,'String',num2str(aux{10}));
    set(handles.value3,'String',num2str(aux{11}));
    set(handles.value4,'String',num2str(aux{12}));
    set(handles.value5,'String',num2str(aux{13}));
    set(handles.value6,'String',num2str(aux{14}));
    set(handles.genrel,'Value',aux{15});
    set(handles.srp,'Value',aux{16});
    set(handles.ealb,'Value',aux{17});  
    set(handles.mass,'String',num2str(aux{18}));
    set(handles.cross_section,'String',num2str(aux{19}));
    set(handles.cr_srp,'String',num2str(aux{20}));
    set(handles.cr_alb,'String',num2str(aux{21}));
    
    handles_array = [handles.mass,handles.cr_alb,handles.cr_srp,handles.cross_section,...
        handles.edit36,handles.edit37,handles.coord1,handles.coord2,handles.coord3,...
        handles.coord4,handles.coord5,handles.coord6,handles.value1,handles.value2,...
        handles.value3,handles.value4,handles.value5,handles.value6,handles.unit1,...
        handles.unit2,handles.unit3,handles.unit4,handles.unit5,handles.unit6,...
        handles.starttime,handles.frame,handles.ctype];
    set(handles_array,'Enable','off');
    
    set(handles.srp,'UserData',aux);
    set(handles.ealb,'UserData',aux);

else
    aux{1} = 1;
    set(handles.srp,'UserData',aux);
    set(handles.ealb,'UserData',aux);
    for j = 1:6
        set(eval(['handles.coord',num2str(j)]),'Enable','inactive');
        set(eval(['handles.unit',num2str(j)]),'Enable','inactive');
    end
    set(handles.mass,'Enable','off');
    set(handles.cr_alb,'Enable','off','String','0');
    set(handles.cr_srp,'Enable','off','String','0');
    set(handles.edit36,'Enable','off');
    set(handles.edit37,'Enable','off');
    set(handles.cross_section,'Enable','off');
end
    ctype_Callback(handles.ctype, eventdata, handles);

% UIWAIT makes computeOrbitGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = computeOrbitGUI_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function starttime_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = str2double(get(hObject,'String'));
if ~isnan(aux)
    aux = cspice_et2utc(aux+64.184,'C',3);
    set(hObject,'String',aux);
end
% Hints: get(hObject,'String') returns contents of starttime as text
%        str2double(get(hObject,'String')) returns contents of starttime as a double


% --- Executes during object creation, after setting all properties.
function starttime_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to starttime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in frame.
function frame_Callback(hObject, eventdata, handles)
% hObject    handle to frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns frame contents as cell array
%        contents{get(hObject,'Value')} returns selected item from frame


% --- Executes during object creation, after setting all properties.
function frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ctype.
function ctype_Callback(hObject, eventdata, handles)
% hObject    handle to ctype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
aux      = contents{get(hObject,'Value')};
if strcmp(aux,'Cartesian')
    coord_array = {'x';'y';'z';'Vx';'Vy';'Vz'};
    unit_array  = {'km';'km';'km';'km/s';'km/s';'km/s'};
    for j = 1:6
        set(eval(['handles.coord',num2str(j)]),'String',coord_array{j});
        set(eval(['handles.unit',num2str(j)]),'String',unit_array{j});
    end
else
    coord_array = {'SMA';'ECC';'INC';'LAN';'AOP';'MA'};
    unit_array  = {'km';' ';'deg';'deg';'deg';'deg'};
    for j = 1:6
        set(eval(['handles.coord',num2str(j)]),'String',coord_array{j});
        set(eval(['handles.unit',num2str(j)]),'String',unit_array{j});
    end
end
% Hints: contents = cellstr(get(hObject,'String')) returns ctype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ctype


% --- Executes during object creation, after setting all properties.
function ctype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ctype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function coord1_Callback(hObject, eventdata, handles)
% hObject    handle to coord1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord1 as text
%        str2double(get(hObject,'String')) returns contents of coord1 as a double


% --- Executes during object creation, after setting all properties.
function coord1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coord1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function value1_Callback(hObject, eventdata, handles)
% hObject    handle to value1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value1 as text
%        str2double(get(hObject,'String')) returns contents of value1 as a double


% --- Executes during object creation, after setting all properties.
function value1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit1_Callback(hObject, eventdata, handles)
% hObject    handle to unit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit1 as text
%        str2double(get(hObject,'String')) returns contents of unit1 as a double


% --- Executes during object creation, after setting all properties.
function unit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sun.
function sun_Callback(hObject, eventdata, handles)
% hObject    handle to sun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sun


% --- Executes on button press in earth.
function earth_Callback(hObject, eventdata, handles)
% hObject    handle to earth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of earth



function degree_Callback(hObject, eventdata, handles)
% hObject    handle to degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
order_Callback(handles.order, eventdata, handles);
% Hints: get(hObject,'String') returns contents of degree as text
%        str2double(get(hObject,'String')) returns contents of degree as a double



function order_Callback(hObject, eventdata, handles)
% hObject    handle to order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = min([str2double(get(handles.degree,'String')),str2double(get(hObject,'String'))]);
if str2double(get(handles.degree,'String'))==1
    set(hObject,'String',num2str(0));
else
    set(hObject,'String',num2str(aux));
end
% Hints: get(hObject,'String') returns contents of order as text
%        str2double(get(hObject,'String')) returns contents of order as a double


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function coord2_Callback(hObject, eventdata, handles)
% hObject    handle to coord2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord2 as text
%        str2double(get(hObject,'String')) returns contents of coord2 as a double



function value2_Callback(hObject, eventdata, handles)
% hObject    handle to value2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value2 as text
%        str2double(get(hObject,'String')) returns contents of value2 as a double


% --- Executes during object creation, after setting all properties.
function value2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit2_Callback(hObject, eventdata, handles)
% hObject    handle to unit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit2 as text
%        str2double(get(hObject,'String')) returns contents of unit2 as a double



function coord3_Callback(hObject, eventdata, handles)
% hObject    handle to coord3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord3 as text
%        str2double(get(hObject,'String')) returns contents of coord3 as a double



function value3_Callback(hObject, eventdata, handles)
% hObject    handle to value3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value3 as text
%        str2double(get(hObject,'String')) returns contents of value3 as a double


% --- Executes during object creation, after setting all properties.
function value3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit3_Callback(hObject, eventdata, handles)
% hObject    handle to unit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit3 as text
%        str2double(get(hObject,'String')) returns contents of unit3 as a double



function coord4_Callback(hObject, eventdata, handles)
% hObject    handle to coord4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord4 as text
%        str2double(get(hObject,'String')) returns contents of coord4 as a double



function value4_Callback(hObject, eventdata, handles)
% hObject    handle to value4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value4 as text
%        str2double(get(hObject,'String')) returns contents of value4 as a double


% --- Executes during object creation, after setting all properties.
function value4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit4_Callback(hObject, eventdata, handles)
% hObject    handle to unit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit4 as text
%        str2double(get(hObject,'String')) returns contents of unit4 as a double



function coord5_Callback(hObject, eventdata, handles)
% hObject    handle to coord5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord5 as text
%        str2double(get(hObject,'String')) returns contents of coord5 as a double



function value5_Callback(hObject, eventdata, handles)
% hObject    handle to value5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value5 as text
%        str2double(get(hObject,'String')) returns contents of value5 as a double


% --- Executes during object creation, after setting all properties.
function value5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit5_Callback(hObject, eventdata, handles)
% hObject    handle to unit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit5 as text
%        str2double(get(hObject,'String')) returns contents of unit5 as a double



function coord6_Callback(hObject, eventdata, handles)
% hObject    handle to coord6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of coord6 as text
%        str2double(get(hObject,'String')) returns contents of coord6 as a double



function value6_Callback(hObject, eventdata, handles)
% hObject    handle to value6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value6 as text
%        str2double(get(hObject,'String')) returns contents of value6 as a double


% --- Executes during object creation, after setting all properties.
function value6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function unit6_Callback(hObject, eventdata, handles)
% hObject    handle to unit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of unit6 as text
%        str2double(get(hObject,'String')) returns contents of unit6 as a double


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux{1}   = str2double(get(handles.degree,'String'));
aux{2}   = str2double(get(handles.order,'String'));
aux{3}   = get(handles.sun,'Value');
aux{4}   = get(handles.earth,'Value');
aux{5}   = get(handles.starttime,'String');
contents = cellstr(get(handles.frame,'String'));
aux{6}   = contents{get(handles.frame,'Value')};
contents = cellstr(get(handles.ctype,'String'));
aux{7}   = contents{get(handles.ctype,'Value')};
aux{8}   = str2double(get(handles.value1,'String'));
aux{9}   = str2double(get(handles.value2,'String'));
aux{10}  = str2double(get(handles.value3,'String'));
aux{11}  = str2double(get(handles.value4,'String'));
aux{12}  = str2double(get(handles.value5,'String'));
aux{13}  = str2double(get(handles.value6,'String'));
aux{14}  = get(handles.genrel,'Value');
aux{15}  = get(handles.srp,'Value');
aux{16}  = get(handles.ealb,'Value');
aux{17}  = str2double(get(handles.mass,'String'));
aux{18}  = str2double(get(handles.cross_section,'String'));
aux{19}  = str2double(get(handles.cr_srp,'String'));
aux{20}  = str2double(get(handles.cr_alb,'String'));
aux{21}  = get(handles.stoptime,'String');
aux{22}  = str2double(get(handles.steptime,'String'));
assignin('base','ORBguidata',aux);
close(gcf);


% --- Executes on button press in genrel.
function genrel_Callback(hObject, eventdata, handles)
% hObject    handle to genrel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of genrel


% --- Executes on button press in srp.
function srp_Callback(hObject, eventdata, handles)
% hObject    handle to srp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = get(hObject,'Value');
if hObject.UserData{1}
    if aux
        set(handles.cr_srp,'Enable','on','String','1.0');
        set(handles.mass,'Enable','on');
        set(handles.cross_section,'Enable','on');    
        set(handles.edit36,'Enable','on');
        set(handles.edit37,'Enable','on');
    else
        set(handles.cr_srp,'Enable','off','String','0');
        aux = get(handles.ealb,'Value');
        if aux
            set(handles.mass,'Enable','on');
            set(handles.cross_section,'Enable','on');
            set(handles.cr_alb,'Enable','on');
            set(handles.edit36,'Enable','on');
            set(handles.edit37,'Enable','on');
        else
            set(handles.mass,'Enable','off');
            set(handles.cross_section,'Enable','off');
            set(handles.edit36,'Enable','off');
            set(handles.edit37,'Enable','off');
            set(handles.cr_alb,'Enable','off','String','0');
        end        
    end
else
    set(handles.mass,'Enable','off','String',hObject.UserData{18});
    set(handles.cross_section,'Enable','off','String',hObject.UserData{19});  
    if aux
        set(handles.cr_srp,'Enable','off','String',hObject.UserData{20});
    else
        set(handles.cr_srp,'Enable','off','String','0');
    end        
    set(handles.edit36,'Enable','off');
    set(handles.edit37,'Enable','off');
end
    
% Hint: get(hObject,'Value') returns toggle state of srp


% --- Executes on button press in ealb.
function ealb_Callback(hObject, eventdata, handles)
% hObject    handle to ealb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = get(hObject,'Value');
if hObject.UserData{1}
    if aux
        set(handles.cr_alb,'Enable','on','String','0.3');
        set(handles.mass,'Enable','on');
        set(handles.cross_section,'Enable','on');  
        set(handles.edit36,'Enable','on');
        set(handles.edit37,'Enable','on');
    else
        set(handles.cr_alb,'Enable','off','String','0');
        aux = get(handles.srp,'Value');
        if aux
            set(handles.mass,'Enable','on');
            set(handles.cross_section,'Enable','on');
            set(handles.edit36,'Enable','on');
            set(handles.edit37,'Enable','on');        
            set(handles.cr_srp,'Enable','on');
        else
            set(handles.mass,'Enable','off');
            set(handles.cross_section,'Enable','off');
            set(handles.edit36,'Enable','off');
            set(handles.edit37,'Enable','off');
            set(handles.cr_srp,'Enable','off','String','0');
        end        
    end
else
    set(handles.mass,'Enable','off','String',hObject.UserData{18});
    set(handles.cross_section,'Enable','off','String',hObject.UserData{19});  
    if aux
        set(handles.cr_alb,'Enable','off','String',hObject.UserData{21});
    else
        set(handles.cr_alb,'Enable','off','String','0');
    end        
    set(handles.edit36,'Enable','off');
    set(handles.edit37,'Enable','off');
end
    
% Hint: get(hObject,'Value') returns toggle state of ealb



function mass_Callback(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass as text
%        str2double(get(hObject,'String')) returns contents of mass as a double


% --- Executes during object creation, after setting all properties.
function mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cross_section_Callback(hObject, eventdata, handles)
% hObject    handle to cross_section (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cross_section as text
%        str2double(get(hObject,'String')) returns contents of cross_section as a double


% --- Executes during object creation, after setting all properties.
function cross_section_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cross_section (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cr_srp_Callback(hObject, eventdata, handles)
% hObject    handle to cr_srp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = str2double(get(hObject,'String'));
if aux>2
    set(hObject,'String','2.0');
end
% Hints: get(hObject,'String') returns contents of cr_srp as text
%        str2double(get(hObject,'String')) returns contents of cr_srp as a double

% --- Executes during object creation, after setting all properties.
function cr_srp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr_srp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cr_alb_Callback(hObject, eventdata, handles)
% hObject    handle to cr_alb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = str2double(get(hObject,'String'));
if aux>1
    set(hObject,'String','1.0');
end

% Hints: get(hObject,'String') returns contents of cr_alb as text
%        str2double(get(hObject,'String')) returns contents of cr_alb as a double


% --- Executes during object creation, after setting all properties.
function cr_alb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr_alb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double


% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stoptime_Callback(hObject, eventdata, handles)
% hObject    handle to stoptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
id = 0;
aux = get(hObject,'String');
if ~isempty(strfind(aux,'days'))&&~isempty(strfind(aux,'+')) %#ok<STREMP>
    aux = strsplit(aux,' ');
    aux = str2double(aux{1})*86400;
    id  = 1;
elseif strfind(aux,'+')>0
    aux = str2double(aux);
    id = 2;
else
    aux = str2double(get(hObject,'String'));
    if isnan(aux)
        id = 3;
    end
end
if id==1||id==2
    aux = cspice_et2utc(aux+cspice_str2et(handles.starttime.String),'C',3);
    set(hObject,'String',aux);
elseif id==0
    aux = cspice_et2utc(aux+64.184,'C',3);
    set(hObject,'String',aux);
end

% Hints: get(hObject,'String') returns contents of stoptime as text
%        str2double(get(hObject,'String')) returns contents of stoptime as a double


% --- Executes during object creation, after setting all properties.
function stoptime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function steptime_Callback(hObject, eventdata, handles)
% hObject    handle to steptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
aux = str2double(get(hObject,'String'));
if isnan(aux)
    set(hObject,'String','60');
end

% Hints: get(hObject,'String') returns contents of steptime as text
%        str2double(get(hObject,'String')) returns contents of steptime as a double


% --- Executes during object creation, after setting all properties.
function steptime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to steptime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
