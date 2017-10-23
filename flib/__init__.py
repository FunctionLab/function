import importlib
import os

ENVIRONMENT_VARIABLE = 'FLIB_SETTINGS'

class Settings:
    def __init__(self, settings_module):

        self.SETTINGS_MODULE = settings_module

        mod = importlib.import_module(self.SETTINGS_MODULE)

        for setting in dir(mod):
            if setting.isupper():
                setting_value = getattr(mod, setting)
                setattr(self, setting, setting_value)

flib_settings = os.environ.get(ENVIRONMENT_VARIABLE, 'flib.settings_default')
settings = Settings(flib_settings)
