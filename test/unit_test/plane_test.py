"""
sample unit test - test travis setup
"""
import sys
import unittest

sys.path.insert(1, '../../src')
import z_plane

class TestPlane(unittest.TestCase):
    def setUp(self):
        self.frame_dict = {'theta': 0.0,
                           'center_point': complex(0.0, 0.0),
                           'zoom_factor': 1,
                           'n_rows': 5,
                           'n_cols': 5 }
        self.frame = {'center_point': complex(0.0, 0.0),
                      'top_center': complex(0.0, 1.0),
                      'right_center': complex(1.0, 0.0),
                      'bottom_center': complex(0.0, -1.0),
                      'left_center': complex(-1.0, 0.0),
                      'upper_right': complex(1.0, 1.0),
                      'bottom_right': complex(1.0, -1.0),
                      'upper_left': complex(-1.0, 1.0),
                      'bottom_left': complex(-1.0, -1.0) }

    def tearDown(self):
        del self.frame_dict
        del self.frame

    def test_get_complex_frame(self):
        CP = self.frame_dict['center_point']
        ZM = self.frame_dict['zoom_factor']
        theta = self.frame_dict['theta']
        frm_def = z_plane.get_complex_frame(CP,ZM,theta)
        self.assertEqual(self.frame['center_point'], frm_def['center_point'])
        self.assertAlmostEqual(self.frame['top_center'], frm_def['top_center'])
        self.assertAlmostEqual(self.frame['right_center'], frm_def['right_center'])
        self.assertAlmostEqual(self.frame['bottom_center'], frm_def['bottom_center'])
        self.assertAlmostEqual(self.frame['left_center'], frm_def['left_center'])
        self.assertAlmostEqual(self.frame['upper_right'], frm_def['upper_right'])
        self.assertAlmostEqual(self.frame['upper_left'], frm_def['upper_left'])
        self.assertAlmostEqual(self.frame['bottom_left'], frm_def['bottom_left'])

if __name__ == '__main__':
    unittest.main()