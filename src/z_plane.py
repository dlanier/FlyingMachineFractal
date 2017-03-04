"""
complex plane functions and class for getting pixels in three different measuring systems


        function framePoints = getFrame(h, w, theta, ZM, CP)
            %Usage: framePoints = obj.getFramePoints;
            %struct of vectors that define the pixels on the complex plane

            % scale and rotate two vectors to the object "zoom" factor
            if w >= h
                framePoints.top_center = exp(1i * (pi / 2 + theta)) / ZM;
                framePoints.right = (w/h) * exp(1i * theta) / ZM;
            else
                framePoints.top_center = (h/w) * exp(1i *(pi/2 + theta))/ZM;
                framePoints.right = exp(1i * theta) / ZM;
            end

            % vector addition == other frame point vectors
            framePoints.bottom_center = -1 * framePoints.top_center;
            framePoints.left = -1 * framePoints.right;
            framePoints.upper_right = framePoints.right +...
                                        framePoints.top_center;
            framePoints.lower_right = framePoints.right +...
                                        framePoints.bottom_center;
            framePoints.upper_left = framePoints.left +...
                                        framePoints.top_center;
            framePoints.lower_left = framePoints.bottom_center -...
                                        framePoints.right;

            % (struct) framePoints: shift all to start at the center point
            nms = fieldnames(framePoints);
            for k = 1:numel(nms)
                framePoints.(nms{k}) = framePoints.(nms{k}) + CP;
            end
        end


    upper_left,      top_center,      upper_right

    left_center,     center_point,    right_center

    bottom_left,     bottom_center,   bottom_right
"""
import numpy as np

print_order = ['upper_left', 'top_center', 'upper_right',
               'left_center', 'center_point', 'right_center',
              'bottom_left', 'bottom_center', 'bottom_right']


def get_complex_frame_string(fd):
    """ get a formatted list of strings """
    frame_string = ''
    row = 0
    for k in print_order:
        r = np.real(fd[k])
        i = np.imag(fd[k])
        row += 1
        if np.mod(row,3) == 0:
            frame_string = frame_string + '\t %s \t %0.9f %0.9f \n'%(k,r,i)
        elif np.mod(row,3) == 1:
            frame_string = frame_string + '%s \t %0.9f %0.9f' % (k, r, i)
        else:
            frame_string = frame_string + '\t %s \t %0.9f %0.9f' % (k, r, i)

    return frame_string

def get_complex_frame(CP, ZM, theta, h=1, w=1):
    """ get the complex numbers at ends and centers of a frame  """
    frame_dict = {'center_point':CP}
    if w >= h:
        frame_dict['top_center'] = np.exp(1j*(np.pi/2 + theta))/ZM
        frame_dict['right_center'] = (w/h) * np.exp(1j * theta) / ZM
    else:
        frame_dict['top_center'] = (h/w) * np.exp(1j*(np.pi/2 + theta)) / ZM
        frame_dict['right_center'] = np.exp(1j * theta) / ZM

    frame_dict['bottom_center'] = frame_dict['top_center'] * -1
    frame_dict['left_center'] = frame_dict['right_center'] * -1
    frame_dict['upper_right'] = frame_dict['right_center'] + frame_dict['top_center']
    frame_dict['bottom_right'] = frame_dict['right_center'] + frame_dict['bottom_center']

    frame_dict['upper_left'] = frame_dict['left_center'] + frame_dict['top_center']
    frame_dict['bottom_left'] = frame_dict['left_center'] + frame_dict['bottom_center']

    return frame_dict

























































""" wolf wolf """