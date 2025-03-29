class DelegatedProperty:
    def __init__(self):
        self.getter = None
        self.setter = None
        self.default = None
        self._doc = None

    def __set_name__(self, owner, name):
        self.name = name
        # Optionally, append a documentation line to the class's __doc__.
        doc_line = f"\n{name}: {self._doc}"
        owner.__doc__ = (owner.__doc__ or "") + doc_line

    def __get__(self, instance, owner):
        if instance is None:
            return self
        try:
            return self.getter(instance)
        except Exception as e:
            if self.default is not None:
                return self.default
            raise AttributeError(f"{owner.__name__} cannot get {self.name}: {e}")

    def __set__(self, instance, value):
        if self.allowed_type is not None:
            # Check allowed type before calling setter.
            if not isinstance(value, self.allowed_type):
                expected = (self.allowed_type.__name__ if isinstance(self.allowed_type, type)
                            else ", ".join(t.__name__ for t in self.allowed_type))
                raise TypeError(f"{type(instance).__name__}.{self.name} only allows type(s) {expected}; "
                                f"got {type(value).__name__} instead.")
        if self.setter is None:
            raise AttributeError(f"{type(instance).__name__} cannot set {self.name} (setter not defined)")
        try:
            self.setter(instance, value)
        except Exception as e:
            raise AttributeError(f"{type(instance).__name__} cannot set {self.name}: {e}")

    def delegates(self, getter, setter=None, default=None, doc=None, allowed_type=None):
        """
        Configure the getter, setter, default, and optional documentation.

        Parameters:
        -----------
        getter: callable
            A function to get the value.
        setter: callable, optional
            A function to set the value.
        default: any, optional
            A default value to return if the getter fails.
        doc: str, optional
            A documentation string for the property.
            If not provided, a default docstring is created based on the default value.
        allowed_type: a type or a tuple of types that are allowed when setting.

        Returns:
        --------
        self: DelegatedProperty
            The current instance for method chaining.

        Example:
        --------
        resolution = DelegatedProperty().delegates(
            getter=lambda self: (bpy.context.scene.render.resolution_x,
                                bpy.context.scene.render.resolution_y),
            setter=lambda self, value: (
                setattr(bpy.context.scene.render, 'resolution_x', value[0]),
                setattr(bpy.context.scene.render, 'resolution_y', value[1])
            ),
            default=(1920, 1080),
            doc="Resolution of the scene in pixels.",
            allowed_type=(tuple, list)
        )
        """
        self.getter = getter
        self.setter = setter
        self.default = default
        self.allowed_type = allowed_type
        allowed_info = f" Allowed type: {allowed_type.__name__}" if allowed_type and isinstance(allowed_type, type) \
                       else f" Allowed types: {[t.__name__ for t in allowed_type]}" if allowed_type else ""
        self._doc = doc if doc is not None else f"Delegated property with default {default}.{allowed_info}"
        return self
    
    def temporary_set(self, instance, value):
        """
        Temporarily set the property value for the given instance.

        Parameters:
        -----------
        instance: object
            The instance to set the property on.
        value: any
            The value to set.

        Returns:
        --------
        None
        """
        if self.setter is None:
            raise AttributeError(f"{type(instance).__name__} cannot set {self.name} (setter not defined)")
        try:
            self.setter(instance, value)
        except Exception as e:
            raise AttributeError(f"{type(instance).__name__} cannot set {self.name}: {e}")